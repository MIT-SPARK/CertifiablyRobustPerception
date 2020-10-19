% A running example for robust single rotation averaging
% Heng Yang
% 10/18/2020

clc
clear
close all
addpath(genpath('../lib'))
addpath(genpath('../solver'))
addpath(genpath('../../mosek'))
addpath(genpath('../../YALMIP')) % use YALMIP for SDP relaxation

%% generate random single rotation averaging problem 
problem.N = 10;
problem.outlierRatio = 0.6;
problem.noiseSigma = 3; % in degrees
problem = gen_single_rotation_averaging(problem);

%% TLS-denseSOS (slow)
% solution = tls_single_rotation_averaging_denseSOS(problem);
% R_err = getAngularError(problem.R_gt, solution.R_est);
% cprintf('keywords','%s using %s: R_err = %g[deg].\n',...
%     problem.type, solution.type, R_err);
    
%% TLS-sparseSOS
solution = tls_single_rotation_averaging_sparseSOS(problem);
R_err = getAngularError(problem.R_gt, solution.R_est);
cprintf('keywords','%s using %s: R_err = %g[deg].\n',...
    problem.type, solution.type, R_err);