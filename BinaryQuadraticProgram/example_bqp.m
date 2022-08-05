%% Example: Dense SDP relaxation for binary quadratic programming (BQP)

clc; clear; close all; restoredefaultpath; % start clean

mosekpath = '../../mosek';
pgdpath   = '../STRIDE';
sdpnalpath  = '../../SDPNAL+v1.0';
utilspath   = '../utils';
addpath(genpath(utilspath))
addpath(genpath('./solvers'))
addpath(genpath('../spotless')) % Use spotless for defining polynomials
addpath('../SDPRelaxations') % implementations for SDP relaxation

%% Generate random binary quadratic program
d       = 20; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
Q       = randn(d,d); Q = Q + Q'; % a random symmetric matrix
c       = randn(d,1);
f       = x'*Q*x + c'*x; % objective function of the BQP
h       = x.^2 - 1; % equality constraints of the BQP (binary variables)
% g       = [x(1)]; % ask the first variable to be positive % if you add
% inequality constraints then you need to change the local search method

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
% problem.inequality      = g;
kappa                   = 2; % relaxation order
[SDP,info]              = dense_sdp_relax(problem,kappa);
SDP.M       = length(info.v); % upper bound on the trace of the moment matrix
% need the following for fast computation in the local search method
info.v      = msspoly2degcoeff(info.v);
info.f      = msspoly2degcoeff(info.f);
info.J      = msspoly2degcoeff(info.J);

%{
%% Solve using MOSEK, should be slow
prob       = convert_sedumi2mosek(SDP.sedumi.At,...
                                  SDP.sedumi.b,...
                                  SDP.sedumi.c,...
                                  SDP.sedumi.K);
addpath(genpath(mosekpath))
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
rmpath(genpath(mosekpath))

figure; bar(eig(Xopt{1})); % if rank = 1, then relaxation is exact/tight
%}


%% solve using stride
addpath(genpath(pgdpath))

pgdopts.pgdStepSize     = 10;
pgdopts.SDPNALpath      = sdpnalpath;
pgdopts.tolADMM         = 10e-5;
pgdopts.phase1          = 1;
pgdopts.rrOpt           = 1:3;
pgdopts.rrFunName       = 'local_search_bqp'; % see solvers/local_search_bqp.m for implementation of local search
pgdopts.rrPar           = info; % need the original POP formulation for local search
pgdopts.maxiterLBFGS    = 1000;
pgdopts.maxiterSGS      = 300;
pgdopts.tolLBFGS        = 1e-12;
pgdopts.tolPGD          = 1e-8;

[outPGD,Xopt,yopt,Sopt]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, [], pgdopts);
time_pgd                    = outPGD.totaltime;
% round solutions and check optimality certificate
res = get_performance_bqp(Xopt,yopt,Sopt,SDP,info,pgdpath);


%% helper functions
function s = msspoly2degcoeff(f)
[~,degmat,coeff,~] = decomp(f);
s.degmat = degmat';
s.coefficient = coeff;
end




