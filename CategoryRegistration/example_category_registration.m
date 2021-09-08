%% Example: Random Category Registration problem and its semidefinite relaxations
%% Heng Yang, July 06, 2021

clc; clear; close all; restoredefaultpath

%% paths to dependencies
spotpath    = '../spotless';
mosekpath   = '../../mosek';
stridepath  = '../STRIDE';
manoptpath  = '../manopt';
sdpnalpath  = '../../SDPNAL+v1.0';
path.stridepath = stridepath;
path.mosekpath  = mosekpath;
path.manoptpath = manoptpath;
addpath('../utils')
addpath('./solvers')

%% choose if run GNC for STRIDE
rungnc      = true;

%% generate random point cloud registration problem
datatype                 = 'car';
problem.outlierRatio     = 0.1;
problem.noiseSigma       = 0.01;
problem.translationBound = 10.0;
switch datatype
    case 'random'
        %%%%%%%%%%%%%% RANDOM data %%%%%%%%%%%%%%
        problem.N                = 10;
        problem.K                = 3;
        problem                  = gen_category_registration(problem);
    case 'car'
        %%%%%%%%%%%%%% PASCAL Car %%%%%%%%%%%%%%
        problem.path             = './data/car.mat';
        problem                  = gen_catreg_pascal_car(problem);
end
lambda                   = 0.5;

%% generate SDP relaxations
addpath(genpath(spotpath))
SDP     = relax_category_registration_v2(problem,'lambda',lambda,'checkMonomials',false);
rmpath(genpath(spotpath))

%% Solve using STRIDE
% Primal initialization
if rungnc
    solution        = gnc_category_registration(problem,SDP,path,'lambda',lambda);
    [R_est,t_est]   = invert_transformation(solution.R_est,solution.t_est);

    v0       = lift_catreg_v2(R_est(:),...
                              t_est(:),...
                              solution.c_est(:),...
                              solution.theta_est(:),...
                              problem.cBound,...
                              problem.translationBound);
    X0          = rank_one_lift(v0);
    gnc.R_err   = getAngularError(problem.R_gt,solution.R_est);
    gnc.t_err   = getTranslationError(problem.t_gt,solution.t_est);
    gnc.c_err   = getTranslationError(problem.c_gt,solution.c_est);
    gnc.time    = solution.time;
    gnc.info    = solution;

else
    X0       = [];
end

% dual initialization
addpath(genpath(spotpath))
chordalSDP     = chordal_relax_category_registration_v2(problem,'lambda',lambda);
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
[Xchordal,ychordal,Schordal,~] = recover_mosek_sol_blk(res,chordalSDP.blk);
S_assm             = catreg_dual_from_chordal_dual_v2(Schordal,problem.N,problem.K);


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
pgdopts.rrFunName       = 'local_search_catreg_v2';
rrPar.blk = SDP.blk; rrPar.translationBound = problem.translationBound;
rrPar.N = problem.N; rrPar.K = problem.K; rrPar.cBound = problem.cBound;
pgdopts.rrPar           = rrPar;
pgdopts.maxiterLBFGS    = 1000;
pgdopts.maxiterSGS      = 1000;
pgdopts.lbfgseps        = false;
pgdopts.S0              = S_assm;


[outPGD,Xopt,yopt,Sopt]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, X0, pgdopts);
rmpath(genpath(manoptpath))

infostride                  = get_performance_catreg_v2(Xopt,yopt,Sopt,SDP,problem,stridepath);
infostride.totaltime        = outPGD.totaltime + time_dualInit;
infostride.time             = [outPGD.totaltime,time_dualInit];
if rungnc 
    infostride.gnc  = gnc; 
    infostride.totaltime = infostride.totaltime + gnc.time;
    infostride.time      = [gnc.time,infostride.time];
end

fprintf('\n\n\n\n\n')
