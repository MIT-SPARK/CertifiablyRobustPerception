% A running example for robust point cloud registration
% Heng Yang
% 10/18/2020

clc
clear
close all
addpath(genpath('../lib'))
addpath(genpath('../solver'))
addpath(genpath('../../YALMIP')) % use YALMIP for SDP relaxation

%% generate random point cloud registration problem
% w.l.o.g. assume the translation is bounded by 1
problem.N = 10;
problem.outlierRatio = 0.5;
problem.translationBound = 1.0;
problem.noiseSigma = 1e-2;
problem = gen_point_cloud_registration(problem);

%% solve TLS robust registration using denseSOS relaxation (slow)
% solution = tls_point_cloud_registration_denseSOS(problem);
% R_err = getAngularError(problem.R_gt, solution.R_est);
% t_err = getTranslationError(problem.t_gt, solution.t_est);
% cprintf('keywords','%s using %s: R_err = %g[deg], t_err = %g.\n',...
%     problem.type, solution.type, R_err, t_err);

%% solve TLS robust registration using sparseSOS relaxation
solution = tls_point_cloud_registration_sparseSOS(problem);
R_err = getAngularError(problem.R_gt, solution.R_est);
t_err = getTranslationError(problem.t_gt, solution.t_est);
cprintf('keywords','%s using %s: R_err = %g[deg], t_err = %g.\n',...
    problem.type, solution.type, R_err, t_err);