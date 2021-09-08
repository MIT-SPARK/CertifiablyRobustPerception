%% Example: Single Rotation Averaging and its semidefinite relaxations
%% Heng Yang, June 29, 2021.

clc; clear; close all; restoredefaultpath

%% paths to dependencies
spotpath    = '../spotless';
stridepath  = '../STRIDE';
manoptpath  = '../manopt';
mosekpath   = '../../mosek'; % for dual initialization of STRIDE
sdpnalpath  = '../../SDPNAL+v1.0'; % for ADMM+

%% choose which solvers to run
rungnc      = true;
addpath('../utils')
addpath('./solvers')

%% generate random single averaging problem
problem.N               = 20;
problem.outlierRatio    = 0.1;
problem.noiseSigma      = 5; % in degrees
problem                 = gen_single_rotation_averaging(problem);

%% generate sparse semidefinite relaxation
addpath(genpath(spotpath))
SDP       = relax_single_rotation_averaging(problem);
fprintf('\n\n\n\n\n')
rmpath(genpath(spotpath))

%% Solve using STRIDE
% Primal Initialization using GNC
if rungnc
    solution    = gnc_single_rotation_averaging(problem);
    v           = lift_sra(solution.R_est(:),solution.theta_est);
    X0          = rank_one_lift(v);
    gnc.time    = solution.time_gnc;
    gnc.R_err   = getAngularError(problem.R_gt,solution.R_est);
    gnc.info    = solution;
else
    X0          = [];
end

% Dual initialization using Chordal relaxation
addpath(genpath(spotpath))
chordalSDP      = chordal_relax_single_rotation_averaging(problem);
rmpath(genpath(spotpath))
prob = convert_sedumi2mosek(chordalSDP.sedumi.At,...
                            chordalSDP.sedumi.b,...
                            chordalSDP.sedumi.c,...
                            chordalSDP.sedumi.K);
addpath(genpath(mosekpath))
time0   = tic;
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 20; % set maximum iterations 20
[~,res] = mosekopt('minimize info',prob,param);
time_dualInit = toc(time0);
[~,~,Schordal,~] = recover_mosek_sol_blk(res,chordalSDP.blk);
S_assm           = sra_dual_from_chordal_dual(Schordal);

% STRIDE main algorithm
addpath(genpath(stridepath))
addpath(genpath(manoptpath))

pgdopts.pgdStepSize     = 10;
pgdopts.SDPNALpath      = sdpnalpath;
pgdopts.maxiterPGD      = 5;

pgdopts.tolADMM         = 1e-6;
pgdopts.maxiterADMM     = 1e3;

pgdopts.rrFunName       = 'local_search_sra';
pgdopts.S0              = S_assm;

[outPGD,Xopt,yopt,Sopt]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, X0, pgdopts);
time_pgd                    = outPGD.totaltime;

infostride                = get_performance_sra(Xopt,yopt,Sopt,SDP,problem,stridepath);
infostride.totaltime      = time_pgd + time_dualInit;
infostride.time           = [time_pgd,time_dualInit];
if rungnc 
    infostride.gnc = gnc; 
    infostride.totaltime  = infostride.totaltime + gnc.time;
    infostride.time       = [infostride.time,gnc.time];
end

rmpath(genpath(manoptpath))
fprintf('\n\n\n\n\n')