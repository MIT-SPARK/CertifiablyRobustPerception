clc
clear
close all

addpath(genpath('./lib'))
n = 6;
A = randn(n,n);
A = A + A';

a = svec(A);

A_ = smat(a);