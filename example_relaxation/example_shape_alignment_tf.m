% A running example for robust shape alignment
% Heng Yang
% 10/18/2020

clc
clear
close all
addpath(genpath('../lib'))
addpath(genpath('../solver'))
addpath(genpath('../../mosek'))
addpath(genpath('../../YALMIP')) % use YALMIP for SDP relaxation

%% generate random shape alignment problem
% w.l.o.g. assume the translation is bounded by 1
problem.N = 10;
problem.centerShape = true;
problem.estimateTranslation = false;
problem.outlierRatio = 0.5;
problem.translationBound = 1;
problem.noiseSigma = 1e-2;
problem = gen_shape_alignment(problem);

%% solve TLS shape alignment with denseSOS
% solution = tls_shape_alignment_tf_denseSOS(problem);
% d_err = abs(problem.s_gt - solution.s_est);
% R_err = getAngularError(problem.R_gt, solution.R_est);
% t_err = getTranslationError(problem.t_gt, solution.t_est);
% cprintf('keywords','%s using %s: s_err = %g, R_err = %g[deg], t_err = %g.\n',...
%     problem.type, solution.type, d_err, R_err, t_err);

%% solve TLS shape alignment with sparseSOS
solution = tls_shape_alignment_tf_sparseSOS(problem);
d_err = abs(problem.s_gt - solution.s_est);
R_err = getAngularError(problem.R_gt, solution.R_est);
t_err = getTranslationError(problem.t_gt, solution.t_est);
cprintf('keywords','%s using %s: s_err = %g, R_err = %g[deg], t_err = %g.\n',...
    problem.type, solution.type, d_err, R_err, t_err);
