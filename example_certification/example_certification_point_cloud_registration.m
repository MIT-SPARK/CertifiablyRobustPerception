% a running example for testing certification algorithm for point cloud registration
% Fast Heuristics: GNC-TLS
% Fast Certification: chordal SOS + Douglas-Rachfold splitting
% Heng Yang
% 10/18/2020

clc
clear
close all
addpath(genpath('../lib'))
addpath(genpath('../solver'))
addpath(genpath('../../YALMIP')) % use YALMIP for SDP relaxation
addpath(genpath('../../SOSTOOLS')) % use the multipoly package in SOSTOOLS to manipulate polynomials

%% generate random point cloud registration problem
% w.l.o.g. assume the translation is bounded by 1
problem.N = 10;
problem.outlierRatio = 0.5;
problem.translationBound = 1.0;
problem.noiseSigma = 0.01;
problem = gen_point_cloud_registration(problem);

%% solve TLS robust registration using GNC
solution = tls_point_cloud_registration_gnc(problem);
R_err = getAngularError(problem.R_gt, solution.R_est);
t_err = getTranslationError(problem.t_gt, solution.t_est);
cprintf('keywords','%s using %s: R_err = %g[deg], t_err = %g.\n',...
    problem.type, solution.type, R_err, t_err);

%% certify the GNC solution
certification = tls_point_cloud_registration_certification_hybrid(problem,solution,...
    'maxIters',5e3,...
    'fixIters',true,...
    'plotSuboptTraj',true);
