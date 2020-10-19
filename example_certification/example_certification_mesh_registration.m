% a running example for testing certification algorithm for mesh
% registration
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
addpath(genpath('../../mosek'))

%% generate random mesh registration problem
problem.N = 10;
problem.outlierRatio = 0.5;
problem.translationBound = 1.0;
problem.pointNoiseSigma = 0.01;
problem.normalNoiseSigma = 0.01; % can make noise larger for normals
problem = gen_mesh_registration(problem);

%% solve using GNC-TLS
solution = tls_mesh_registration_gnc(problem);
R_err = getAngularError(problem.R_gt, solution.R_est);
t_err = getTranslationError(problem.t_gt, solution.t_est);
cprintf('keywords','%s using %s: R_err = %g[deg], t_err = %g.\n',...
    problem.type, solution.type, R_err, t_err);

%% certify GNC-TLS
certification = tls_mesh_registration_certification_hybrid(problem,solution,...
    'maxIters',5e3,...
    'fixIters',true,...
    'plotSuboptTraj',true);

