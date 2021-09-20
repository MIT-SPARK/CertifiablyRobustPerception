%% Example: Dense SDP relaxation for binary quadratic programming (BQP)

clc; clear; close all; restoredefaultpath; % start clean

mosekpath = '../mosek';

addpath(genpath('./utils'))
addpath(genpath('./spotless')) % Use spotless for defining polynomials
addpath('./SDPRelaxations') % implementations for SDP relaxation techniques

%% Generate random binary quadratic program
d       = 10; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
Q       = randn(d,d); Q = Q + Q'; % a random symmetric matrix
c       = randn(d,1);
f       = x'*Q*x + c'*x; % objective function of the BQP
h       = x.^2 - 1; % equality constraints of the BQP (binary variables)
g       = [x(1)]; % ask the first variable to be positive

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
kappa                   = 2; % relaxation order
[SDP,info]              = dense_sdp_relax(problem,kappa);

%% Solve using MOSEK
prob       = convert_sedumi2mosek(SDP.sedumi.At,...
                                  SDP.sedumi.b,...
                                  SDP.sedumi.c,...
                                  SDP.sedumi.K);
addpath(genpath(mosekpath))
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
rmpath(genpath(mosekpath))

figure; bar(eig(Xopt{1})); % if rank = 1, then relaxation is exact/tight

