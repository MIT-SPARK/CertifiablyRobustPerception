% a running example for testing certification algorithm for single rotation averaging
% Fast Heuristics: GNC-TLS
% Fast Certification: chordal SOS + Douglas-Rachfold splitting
% Heng Yang
% 10/18/2020

clc
clear
close all
addpath(genpath('../lib'))
addpath(genpath('../solver'))
addpath(genpath('../../mosek'))
addpath(genpath('../../YALMIP')) % use YALMIP for SDP relaxation
addpath(genpath('../../SOSTOOLS')) % use the multipoly package in SOSTOOLS to manipulate polynomials

%% generate random single rotation averaging problem
problem.N = 10;
problem.outlierRatio = 0.5;
problem.noiseSigma = 5; % in degrees
problem = gen_single_rotation_averaging(problem);

%% solve using GNC-TLS
solution = tls_single_rotation_averaging_gnc(problem);
R_err = getAngularError(problem.R_gt, solution.R_est);
cprintf('keywords','%s using %s: R_err = %g[deg].\n',...
    problem.type, solution.type, R_err);

%% certify GNC-TLS
certification = tls_single_rotation_averaging_certification_hybrid(problem,solution,...
    'maxIters',5e3,...
    'fixIters',true,...
    'plotSuboptTraj',false);

