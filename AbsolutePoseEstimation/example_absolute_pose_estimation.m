%% Example: random Absolute Pose Estimation problem and its semidefinie relaxation
%% Also known as perspective-n-points (PnP) problem
%% Heng Yang, August 05, 2021

clc; clear; close all; restoredefaultpath

%% paths to dependencies
spotpath    = '../spotless';
mosekpath   = '../../mosek';
yalmippath  = '../../YALMIP';
stridepath  = '../STRIDE';
manoptpath  = '../manopt';
sdpnalpath  = '../../SDPNAL+v1.0';
pcrpath     = '../PointCloudRegistration/solvers';
addpath('../utils')
addpath('./solvers')

%% choose whether to run GNC
rungnc      = true;

%% generate random absolute pose estimation problem
problem.N                = 10;
problem.outlierRatio     = 0.1;
problem.noiseSigma       = 0.001;
problem.translationBound = 20.0;
problem                  = gen_absolute_pose_estimation(problem);

%% relax the problem into an SDP
addpath(genpath(spotpath))
SDP       = relax_absolute_pose_estimation_v1(problem,'checkMonomials',false);
fprintf('\n\n\n\n\n')
rmpath(genpath(spotpath))

%% Solving using STRIDE
addpath(genpath(stridepath))
addpath(genpath(manoptpath))
addpath(genpath(pcrpath))
% primal initialization
if rungnc
    solution            = msac_absolute_pose_estimation(problem,SDP);
    gnc.time            = solution.time;
    gnc.f_est           = solution.f_est;
    gnc.info            = solution;
    gnc.R_err           = getAngularError(problem.R_gt,solution.R_est);
    gnc.t_err           = getTranslationError(problem.t_gt,solution.t_est);
    v                   = lift_ape_v1(solution.R_est(:),solution.t_est,solution.theta_est,...
                           problem.translationBound,problem.depthBound,problem.FOV);
    X0                  = rank_one_lift(v);
else
    X0       = [];
end

% dual initialization
addpath(genpath(spotpath))
chordalSDP       = chordal_relax_absolute_pose_estimation(problem);
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
[~,~,Sopt,obj] = recover_mosek_sol_blk(res,chordalSDP.blk);
S_assm               = ape_dual_from_chordal_dual(Sopt);

% start STRIDE
pgdopts.pgdStepSize     = 10;
pgdopts.maxiterPGD      = 5;
pgdopts.SDPNALpath      = sdpnalpath;
pgdopts.tolADMM         = 1e-12;
pgdopts.maxiterADMM     = 1e4;
pgdopts.stopoptionADMM  = 0;

pgdopts.rrFunName       = 'local_search_ape_v1';

rrPar.blk               = SDP.blk; 
rrPar.translationBound  = problem.translationBound;
rrPar.depthBound        = problem.depthBound;
rrPar.FOV               = problem.FOV;

pgdopts.rrPar           = rrPar;
pgdopts.maxiterLBFGS    = 1000;
pgdopts.maxiterSGS      = 1000;
pgdopts.lbfgseps        = false;
pgdopts.S0              = S_assm;

[outPGD,Xopt,yopt,Sopt] = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, X0, pgdopts);

rmpath(genpath(manoptpath))
rmpath(genpath(pcrpath))

infostride              = get_performance_ape_v1(Xopt,yopt,Sopt,SDP,problem,stridepath);
infostride.totaltime    = outPGD.totaltime + time_dualInit;
infostride.time         = [outPGD.totaltime,time_dualInit];
if rungnc
    infostride.gnc          = gnc; 
    infostride.totaltime    = infostride.totaltime + gnc.time;
    infostride.time         = [infostride.time, gnc.time];
end
fprintf('\n\n\n\n\n')
