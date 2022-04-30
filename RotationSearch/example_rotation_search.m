%% Robust rotation search (Wahba Problem)
%% Yang, Heng, and Luca Carlone, ICCV 2019
%% "A quaternion-based certifiably optimal solution to the Wahba problem with outliers." 
%% Heng Yang, April 30, 2022

clc
clear
close all
restoredefaultpath

cvxpath     = '../../cvx';
mosekpath   = '../../mosek';
utilspath   = '../utils';
addpath(genpath(pwd))
addpath(genpath(cvxpath))
addpath(genpath(mosekpath))
addpath(genpath(utilspath))

%% generate random problem
N       = 20; % number of measurements
outrate = 0.5; % outlier rate
problem.N               = N;
problem.Covariance      = 1e-4*eye(3);
problem.v1_distribution = 'uniform';
problem.nrOutliers      = round(N*outrate);
problem.boundNoise      = true;
problem.normalize       = false;

[a, b, R_gt, problem]   = createWahbaProblem(problem);
betasq = 1;
if problem.boundNoise; betasq = betasq*(problem.noiseBound)^2; end
if betasq == 0; betasq = 1e-2; end

%% solve SDP relaxation
C = cost(a,b,betasq);
n = 4*N + 4;
cvx_begin % sdp
cvx_solver mosek
variable X(n,n) symmetric
minimize( trace(C * X) ) %
subject to
X == semidefinite(n);
trace(X((1:4),(1:4))) == 1;
for k=1:N
    idx = 4+blkIndices(k,4);
    X(idx,idx) == X((1:4),(1:4));
end
for k1=1:N
    for k2=k1+1:N+1
        idx1 = blkIndices(k1,4);
        idx2 = blkIndices(k2,4);
        X(idx1,idx2) == X(idx1,idx2)';
    end
end
cvx_end

%% extract solution
f_sdp   = cvx_optval; % lower bound
[V,~]   = sorteig(X);
v       = V(:,1);
q       = normalize_cols( v(1:4) );
theta   = zeros(N,1);
for i = 1:N
    theta(i) = sign(q'*v(blkIndices(i+1,4)));
end
R_err   = getAngularError(quat2rot(q),R_gt);
f_est   = cost_org(a,b,betasq,q); % upper bound
subopt  = abs(f_est - f_sdp) / (1+abs(f_est)+abs(f_sdp));

fprintf('f_sdp: %3.4e, f_est: %3.4e, R_err: %3.2e[deg].\n',f_sdp,f_est,R_err);
fprintf('Relative suboptimality: %3.2e.\n',subopt);


%% *********************************************
%% *********************************************
%% helper function
function Q_cost = cost(v1,v2,barc2)
P=[1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
   0, 0, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0;
   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0;
   -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
   0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0;
   0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
   0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0;
   -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
P=sparse(P);
N = size(v1,2);
Npm = 4*N + 4;
Q_1=zeros(Npm,Npm);
for k=1:N
    idx = 4+blkIndices(k,4);
    P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
    ck = 0.5 * ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) - barc2 );
    Q_1((1:4),idx) = Q_1((1:4),idx)-0.5*P_k+ck/2*eye(4);
    Q_1(idx,(1:4)) = Q_1(idx,(1:4))-0.5*P_k+ck/2*eye(4);
end

Q_2=zeros(Npm,Npm);
for k=1:N
%     idx = 4+blkIndices(k,4);
    idx = blkIndices(1,4);
    P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
    ck = 0.5 * ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) + barc2 );
    Q_2(idx,idx) = Q_2(idx,idx) - P_k + ck*eye(4);
end
Q_cost=Q_1+Q_2;
Q_cost=sparse(Q_cost);

Q_cost=Q_cost/barc2;
end

function B = normalize_cols(A)
mag = sum(A.^2,1).^(0.5);
B   = A./mag;
end

function f_est = cost_org(a,b,betasq,q)
R         = quat2rot(q);
residuals = sum((b - R*a).^2,1)/betasq;
f_est     = sum( min(residuals(:),1) );
end

