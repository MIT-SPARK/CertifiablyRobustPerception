clc
clear
close all
restoredefaultpath;

mosekpath = '../../mosek';
utilspath   = '../utils';
addpath(genpath(utilspath))
addpath(genpath('../spotless')) % Use spotless for defining polynomials
addpath('../SDPRelaxations') % implementations for SDP relaxation
addpath(genpath('../STRIDE'))
addpath(genpath(mosekpath))

N = 10;
Y = randn(3,N);
R = eye(3);
t = [0;0;10];
Yp = R * Y + t;
y  = Yp(1:2,:) ./ Yp(3,:) + 0.02 * randn(2,N);

[R_est,t_est,info] = SolveWeightedPnP(y,Y,eye(3),ones(N,1));

R_err = getAngularError(R,R_est);
t_err = getTranslationError(t,t_est);

fprintf('R_err: %3.2e, t_err: %3.2e.\n',R_err,t_err);

SDP = info.SDP;
At  = SDP.At{1};
m   = size(SDP.b,1);
n   = SDP.blk{1,2};

A = zeros(m,n,n);
for i = 1:m
    A(i,:,:) = smat(SDP.blk,At(:,i));
end


function [R_est,t_est,info] = SolveWeightedPnP(y,Y,K,c)
%% Inputs
% y: 2 x N image pixels
% Y: 3 x N 3D points
% K: camera intrinsics
% c: N x 1 weights
%% Outputs
% R, t: 6D pose
% info: relaxation tightness

%% preprocess data
N = size(y,2);
yh = [y;ones(1,N)];
tv = linsolve(K,yh);
v  = tv ./ (sum(tv.^2,1).^(1/2));
W  = zeros(N,3,3);
cW = zeros(N,3,3);
cWsum = zeros(3,3);
for i = 1:N
    vi = v(:,i);
    Wi = eye(3) - vi * vi';
    cWi = c(i) * Wi;
    
    W(i,:,:) = Wi;
    cW(i,:,:) = cWi;
    cWsum = cWsum + cWi;
end
bW = zeros(N,3,3);
for i = 1:N
    bWi = linsolve(cWsum,squeeze(cW(i,:,:)));
    bW(i,:,:) = bWi;
end
M = zeros(3,9);
for i = 1:N
    bWi = squeeze(bW(i,:,:));
    Yi  = Y(:,i);
    M = M + kron(Yi',bWi);
end
Q = zeros(9,9);
for i = 1:N
    Yi = Y(:,i);
    Gi = kron(Yi',eye(3)) - M;
    cWi = squeeze(cW(i,:,:));
    Q = Q + Gi' * cWi * Gi;
end

%% define POP
d       = 9; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
R       = reshape(x,3,3);
c1      = R(:,1);
c2      = R(:,2);
c3      = R(:,3);
f       = x' * Q * x;
h       = [c1'*c1 - 1;
           c2'*c2 - 1;
           c3'*c3 - 1;
           c1'*c2;
           c2'*c3;
           c3'*c1;
           cross(c1,c2) - c3;
           cross(c2,c3) - c1;
           cross(c3,c1) - c2];
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
kappa                   = 2; % relaxation order
[SDP,relaxinfo]              = dense_sdp_relax(problem,kappa);

prob       = convert_sedumi2mosek(SDP.sedumi.At,...
                                  SDP.sedumi.b,...
                                  SDP.sedumi.c,...
                                  SDP.sedumi.K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);

Xopt = Xopt{1};
[V,~] = sorteig(Xopt);
v1 = V(:,1) / V(1,1);
R_est = project2SO3(reshape(v1(2:10),3,3));
t_est = zeros(3,1);
for i = 1:N
    bWi = squeeze(bW(i,:,:));
    Yi  = Y(:,i);
    t_est = t_est - bWi * R_est * Yi;
end

f_sdp = obj(1);
r_est = R_est(:);
f_est = r_est' * Q * r_est;
eta   = abs(f_est - f_sdp) / (1 + abs(f_est) + abs(f_sdp));

fprintf("f_sdp: %3.2e, f_est: %3.2e, gap: %3.2e.\n",f_sdp,f_est,eta);

info.SDP = SDP;
info.f_est = f_est;
info.f_sdp = f_sdp;
info.eta = eta;
end
