%% Example: mesh registration and its semidefinite relaxation

clc; clear; close all; restoredefaultpath

%% path to dependencies
spotpath    = '../spotless';
mosekpath   = '../../mosek';
yalmippath  = '../../YALMIP';
stridepath  = '../STRIDE';
manoptpath  = '../manopt';
sdpnalpath  = '../../SDPNAL+v1.0';
pcrpath     = '../PointCloudRegistration/solvers';
mypath      = {stridepath,pcrpath};

%% Choose whether to run GNC
rungnc      = true;

addpath('../utils')
addpath('./solvers')

%% generate random point cloud registration problem
problem.N                   = 10;
problem.outlierRatio        = 0.1;
problem.pointNoiseSigma     = 0.01;
problem.normalNoiseSigma    = 0.01;
problem.translationBound    = 10.0;
problem                     = gen_mesh_registration(problem);

%% generate SDP relaxations
addpath(genpath(spotpath))
SDP       = relax_mesh_registration(problem,'checkMonomials',false);
fprintf('\n\n\n\n\n')
rmpath(genpath(spotpath))

%% STRIDE
addpath(genpath(pcrpath))
% primal initialization
if rungnc
    addpath(genpath(yalmippath))
    addpath(genpath(mosekpath))
    solution = gnc_mesh_registration(problem);
    rmpath(genpath(mosekpath))
    rmpath(genpath(yalmippath))
    v        = lift_pcr_v4(solution.R_est(:),...
                           solution.t_est,...
                           solution.theta_est,...
                           problem.translationBound);
    X0       = rank_one_lift(v);

    gnc.R_err = getAngularError(problem.R_gt,solution.R_est_org);
    gnc.t_err = getTranslationError(problem.t_gt,solution.t_est_org);
    gnc.time  = solution.time_gnc;
    gnc.f_est = solution.f_est;
    gnc.info  = solution;
else
    X0       = [];
end

% Dual initialization using chordal SDP
addpath(genpath(spotpath))
chordalSDP       = chordal_relax_mesh_registration(problem);
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
[Xchordal,ychordal,Schordal,~] = recover_mosek_sol_blk(res,chordalSDP.blk);
S_assm           = pcr_dual_from_chordal_dual(Schordal);

% STRIDE main algorithm
addpath(genpath(stridepath))
addpath(genpath(manoptpath))

pgdopts.pgdStepSize     = 10;
pgdopts.maxiterPGD      = 5;
pgdopts.SDPNALpath      = sdpnalpath;
pgdopts.S0              = S_assm;

pgdopts.tolADMM         = 1e-8;
pgdopts.maxiterADMM     = 1e4;
pgdopts.stopoptionADMM  = 0;

pgdopts.rrFunName       = 'local_search_pcr_v4';
rrPar.blk = SDP.blk; rrPar.translationBound = problem.translationBound;
pgdopts.rrPar           = rrPar;
pgdopts.maxiterLBFGS    = 1000;
pgdopts.maxiterSGS      = 1000;
pgdopts.lbfgseps        = false;

[outPGD,Xopt,yopt,Sopt]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, X0, pgdopts);

infostride                  = get_performance_mr(Xopt,yopt,Sopt,SDP,problem,mypath);
infostride.totaltime        = outPGD.totaltime + time_dualInit;
infostride.time             = [outPGD.totaltime,time_dualInit];
if rungnc
    infostride.gnc          = gnc; 
    infostride.totaltime    = infostride.totaltime + gnc.time;
    infostride.time         = [infostride.time, gnc.time];
end

rmpath(genpath(manoptpath))
fprintf('\n\n\n\n\n')