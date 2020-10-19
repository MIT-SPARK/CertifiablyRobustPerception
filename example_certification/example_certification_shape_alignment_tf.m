% a running example for testing certification algorithm for shape alignment
% translation free
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

%% generate random shape alignment problem
problem.N = 10;
problem.outlierRatio = 0.4;
problem.estimateTranslation = false;
problem.noiseSigma = 1e-2;
problem = gen_shape_alignment(problem);


%% solve using GNC-TLS
solution = tls_shape_alignment_tf_gnc(problem);
s_err = abs(problem.s_gt - solution.s_est);
R_err = getAngularError(problem.R_gt, solution.R_est);
t_err = getTranslationError(problem.t_gt, solution.t_est);
cprintf('keywords','%s using %s: s_err = %g, R_err = %g[deg], t_err = %g.\n',...
    problem.type, solution.type, s_err, R_err, t_err);

%% certify the solution
certification = tls_shape_alignment_tf_certification_hybrid(problem,solution,...
    'maxIters',5e3,...
    'fixIters',true,...
    'plotSuboptTraj',true);

