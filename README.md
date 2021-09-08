# Certifiable Outlier-Robust Geometric Perception

## About
This repository holds the implementation for certifiably solving outlier-robust geometric perception problems to global optimality. The certifiable outlier-robust geometric perception framework contains two main modules:

- A sparse semidefinite programming relaxation (SSR) scheme that relaxes nonconvex outlier-robust perception problems into convex semidefinite programs (SDPs); and

- A novel SDP solver, called [STRIDE](https://github.com/MIT-SPARK/STRIDE), that solves the generated SDPs at an unprecedented scale and accuracy.

We proposed a preliminary version of the SSR scheme in our [NeurIPS 2020 paper](https://arxiv.org/abs/2006.06769), and released a certifier (that certifies if a given estimate is optimal) based on Douglas-Rachford Splitting (DRS). Please switch to the `NeurIPS2020` branch of this repo to checkout the NeurIPS2020 implementation.