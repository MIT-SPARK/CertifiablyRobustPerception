%% Example: Random Point Cloud Registration problem and its semidefinite relaxations
%% Heng Yang, June 29, 2021

clc; clear; close all; restoredefaultpath

%% paths to dependencies
spotpath    = '../spotless';
stridepath  = '../STRIDE';
manoptpath  = '../manopt';
mosekpath   = '../../mosek';
sdpnalpath  = '../../SDPNAL+v1.0';
addpath('../utils')
addpath('./solvers')

%% choose if run GNC for STRIDE
rungnc      = true;

%% generate random point cloud registration problem
problem.N                = 10;
problem.outlierRatio     = 0.2;
problem.noiseSigma       = 0.01;
problem.translationBound = 10.0;
problem                  = gen_point_cloud_registration(problem);

%% generate SDP relaxation
addpath(genpath(spotpath))
SDP       = relax_point_cloud_registration_v4(problem,'checkMonomials',false);
fprintf('\n\n\n\n\n')
rmpath(genpath(spotpath))

%% Solve using STRIDE
% primal initialization using GNC
if rungnc
    solution = gnc_point_cloud_registration(problem);
    v        = lift_pcr_v4(solution.R_est(:),...
                           solution.t_est,...
                           solution.theta_est,...
                           problem.translationBound);
    X0       = rank_one_lift(v);

    gnc.R_err = getAngularError(problem.R_gt,solution.R_est);
    gnc.t_err = getTranslationError(problem.t_gt,solution.t_est);
    gnc.time  = solution.time_gnc;
    gnc.f_est = solution.f_est;
    gnc.info  = solution;
else
    X0       = [];
end
% Dual initialization using chordal SDP
addpath(genpath(spotpath))
chordalSDP       = chordal_relax_point_cloud_registration(problem);
fprintf('\n\n\n\n\n')
rmpath(genpath(spotpath))
prob = convert_sedumi2mosek(chordalSDP.sedumi.At,...
                            chordalSDP.sedumi.b,...
                            chordalSDP.sedumi.c,...
                            chordalSDP.sedumi.K);
addpath(genpath(mosekpath))
time0   = tic;
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 20;
[~,res] = mosekopt('minimize info',prob,param);
time_dualInit = toc(time0);
[~,~,Schordal,~] = recover_mosek_sol_blk(res,chordalSDP.blk);
S_assm           = pcr_dual_from_chordal_dual(Schordal);

% STRIDE main algorithm
addpath(genpath(stridepath))
addpath(genpath(manoptpath))

pgdopts.pgdStepSize     = 10;
pgdopts.SDPNALpath      = sdpnalpath;
pgdopts.maxiterPGD      = 5;
% ADMM parameters
pgdopts.tolADMM         = 1e-10;
pgdopts.maxiterADMM     = 1e4;
pgdopts.stopoptionADMM  = 0;

pgdopts.rrOpt           = 1:3;
pgdopts.rrFunName       = 'local_search_pcr_v4';
rrPar.blk = SDP.blk; rrPar.translationBound = problem.translationBound;
pgdopts.rrPar           = rrPar;
pgdopts.maxiterLBFGS    = 1000;
pgdopts.maxiterSGS      = 1000;
pgdopts.S0              = S_assm;

[outPGD,Xopt,yopt,Sopt]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, X0, pgdopts);
rmpath(genpath(manoptpath))

infostride              = get_performance_pcr(Xopt,yopt,Sopt,SDP,problem,stridepath);
infostride.totaltime    = outPGD.totaltime + time_dualInit;
infostride.time         = [outPGD.totaltime,time_dualInit];
if rungnc
    infostride.gnc = gnc; 
    infostride.totaltime = infostride.totaltime + gnc.time;
    infostride.time = [infostride.time, gnc.time];
end
fprintf('\n\n\n\n\n')