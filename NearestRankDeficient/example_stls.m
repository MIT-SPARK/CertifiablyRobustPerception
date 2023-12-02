clc;
clear;
close all;
restoredefaultpath;

sdpnalpath  = '../../SDPNAL+v1.0';
pgdpath     = '../STRIDE';
utilspath   = '../utils';

addpath(genpath(pwd));
addpath(genpath(utilspath))
%% construct space of n1 x n2 hankel matrices (n1 <= n2)
n1 = 20;
n2 = 20;
[S,k,Scell] = hankel_struct(n1,n2);
% generate random hankel matrix
u1 = randn(k,1);
U1 = applyAffineMapCell(Scell,u1);

%% Convert the nearest hankel matrix problem to SDP
addpath(genpath(pgdpath));
SDP     = nearest_hankel_sdp(S,u1);
fprintf('SDP size: n = %d, m = %d.\n\n\n',SDP.n,SDP.m);
rmpath(genpath(pgdpath));


%% Optimistic initialization using SLRA
s.m             = n1;
s.n             = n2;
tol             = 0;
[zhtls,Hu]      = htls(s.m-1,u1,'linear'); % good initialization
opts.Rini       = zhtls';
opts.maxiter    = 5000;
opts.tol        = tol; 
opts.epsrel     = tol;
opts.epsabs     = tol; 
opts.epsgrad    = tol;
opts.method     = 'll'; % 'll': Levenberg–Marquardt, 'qb': BFGS
[uslra, info]   = slra(u1, s, s.m-1, opts);
zslra           = info.Rh'; zslra = zslra/norm(zslra);
Uslra           = applyAffineMapCell(Scell,uslra);
ztUnorm         = norm(zslra'*Uslra); % measure of rank deficientness
fprintf('SLRA: norm(zt*U) = %3.2e.\n',ztUnorm);
% lift to SDP solution
xtld            = kron([uslra;1],zslra);
X0              = {xtld * xtld'};

%% STRIDE
addpath(genpath(pgdpath))
pgdopts.pgdStepSize     = 10;
pgdopts.maxiterSGS      = 300;
pgdopts.maxiterLBFGS    = 1000;
pgdopts.SDPNALpath      = sdpnalpath;
pgdopts.tolADMM         = 5e-5;
pgdopts.phase1          = 1;
pgdopts.rrOpt           = 1:3;
pgdopts.lbfgsmemory     = 10;
pgdopts.tolPGD          = 1e-8;
pgdopts.rrFunName       = 'rr_hankel';
rrPar.m = n1; rrPar.n = n2; rrPar.k = k; rrPar.theta = u1; rrPar.S = S;
pgdopts.rrPar           = rrPar;
[outPGD,Xopt,yopt,Sopt] = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, X0, pgdopts);
f_sdp                   = outPGD.pobj;
time_pgd                = outPGD.totaltime;

[z,u,f_est]     = recover_solution(Xopt,SDP.C,n1,k);
eta             = abs(f_est-f_sdp) / (1+abs(f_est)+abs(f_sdp)); % relative suboptimality gap
U               = applyAffineMapCell(Scell,u);
ztUnorm         = norm(z'*U); % measure of rank deficientness
fprintf('norm(zt*U) = %3.2e, eta = %3.2e.\n',ztUnorm,eta);

rmpath(genpath(pgdpath))
fprintf('\n\n\n\n\n')
