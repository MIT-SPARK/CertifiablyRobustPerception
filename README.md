# Certifiable Outlier-Robust Geometric Perception

## About
This repository holds the implementation for certifiably solving outlier-robust geometric perception problems to global optimality. The certifiable outlier-robust geometric perception framework contains two main modules:

- A sparse semidefinite programming relaxation (SSR) scheme that relaxes nonconvex outlier-robust perception problems into convex semidefinite programs (SDPs); and

- A novel SDP solver, called [STRIDE](https://github.com/MIT-SPARK/STRIDE), that solves the generated SDPs at an unprecedented scale and accuracy.

We proposed a preliminary version of the SSR scheme in our [NeurIPS 2020 paper](https://arxiv.org/abs/2006.06769), and released a certifier (that certifies if a given estimate is optimal) based on Douglas-Rachford Splitting (DRS). Please switch to the `NeurIPS2020` branch of this repo to checkout the NeurIPS2020 implementation.

If you find this library helpful or use it in your projects, please cite:
```bibtex
@article{Yang22tpami-certifiableperception,
  title={Certifiably Optimal Outlier-Robust Geometric Perception: Semidefinite Relaxations and Scalable Global Optimization},
  author={Yang, Heng and Carlone, Luca},
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year={2022}
}
```

## Tutorial: Semidefinite Relaxation for Polynomial Optimization
A general polynomial optimization problem (POP) is an optimization problem of the following standard form
```
minimize_{x ∈ R^d}      f(x)
subject to              h_i(x) == 0, i=1,...,l_h
                        g_j(x) >= 0, j=1,...,l_g
```
where `x ∈ R^d` is a `d`-dimensional decision variable, `f(x)` is a scalar polynomial objective function, `h_i(x),i=1,...,l_h` are scalar polynomial functions that define equality constraints, and `g_j(x),j=1,...,l_g` are polynomial functions that define inequality constraints. POPs are in general nonconvex and NP-hard problems. For example, one can easily see that binary quadratic programming is an instance of POP, because binary constraints `x_i ∈ {+1, -1}, i=1,...,d` can be easily written as polynomial equalities `h_i(x) = x_i^2 - 1 = 0,i=1,...,d`. 

Semidefinite relaxations are a powerful tool for solving nonconvex POPs to **global optimality**. In this repo, we provide tools that can help the user exploit the power of semidefinite relaxations with minimum efforts.

### Dense Relaxation (Lasserre's Hierarchy)

Lasserre's hierarchy of moment relaxations is a standard technique for relaxing POPs into semidefinite programs (SDPs). The basic idea is as follows. Let `kappa` be a positive integer (called the **relaxation order**) such that `2*kappa` is no smaller than the maximum degree of the defining polynomials `f,h_i,g_j`, and let `v = [x]_kappa` be the set of standard monomials in `x` with degree up to `kappa`. For example, suppose `x = [x_1; x_2]`, then `[x]_kappa` with `kappa = 2` leads to `v = [1; x_1; x_2; x_1*x_2; x_1^2; x_2^2]`. The essential idea of Lasserre's hierarchy is to express the original POP problem in the matrix variable `X = v*v'` (called the **moment matrix**, by construction `X` is positive semidefinite) and relax the POP into a convex SDP. In the [seminal paper of Lasserre](https://epubs.siam.org/doi/abs/10.1137/S1052623400366802?journalCode=sjope8), it was proved that if `kappa` is chosen large enough, then the convex SDP can solve the original nonconvex POP to global optimality. 

Although the underlying mechanics of Lasserre's hierarchy can be complicated (the interested reader can refer to Section 2 of our [paper](https://arxiv.org/abs/2109.03349) for technical details), in this repo we provide a simple function that implements Lasserre's hierarchy:
```
[SDP,info] = dense_sdp_relax(problem,kappa)
```
where the inputs are:
- `problem`: a Matlab structure that contains the following fields:
    - `vars`: a vector of symbolic decision variables (i.e., `x` in the POP);
    - `objective`: a multivariate polynomial that specifies the objective function of the POP (i.e., `f(x)` in the POP);
    - `equality` (optional): a vector of multivariate polynomials with dimension `l_h` that specifies the equality constraints of the POP (i.e., `h_i(x),i=1,...,l_h` in the POP);
    - `inequality` (optional): a vector of multivariate polynomials with dimension `l_g` that specifies the inequality constraints of the POP (i.e., `g_j(x),j=1,...,l_g` in the POP).
- `kappa` (optional): a positive integer that specifies the relaxation order. If `kappa` is not provided, or provided such that `2*kappa` is smaller than the maximum degree of the defining polynomials, then our implementation uses the minimum relaxation order `kappa` such that `2*kappa` is no smaller than the maximum degree.

and the outputs are:
- `SDP`: a Matlab structure that contains the following fields:
    - `blk,At,b,C`: standard SDP problem data in SDPT3 format. The interested reader should check out [SDPT3 user manual](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/) for details. Section 2.1 of our [paper](https://arxiv.org/abs/2109.03349) also gives a quick introduction.
    - `sedumi`: a Matlab structure that provides the standard SDP problem data in sedumi format. The interested reader should checkout [sedumi user manual](https://sedumi.ie.lehigh.edu/?page_id=58) for details. It is  easy to convert sedumi format to MOSEK format, as we will show in one of the examples later.
- `info`: a Matlab structure that contains additional information about the relaxation.

#### Example: Binary Quadratic Programming
We now use a simple example on binary quadratic programming (BQP) to illustrate how to supply the POP problem description to `dense_sdp_relax`. The sample code is given below.
```
%% Generate random binary quadratic program
d       = 10; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
Q       = randn(d,d); Q = Q + Q'; % a random symmetric matrix
c       = randn(d,1);
f       = x'*Q*x + c'*x; % objective function of the BQP
h       = x.^2 - 1; % equality constraints of the BQP (binary variables)
g       = [x(1)]; % ask the first variable to be positive

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
kappa                   = 2; % relaxation order
[SDP,info]              = dense_sdp_relax(problem,kappa);
```
In this demo code, we first generate a random binary quadratic programming problem using the package SPOTLESS (which is a submodule of this repo), and then pass the `problem` structure with fields `vars`, `objective`, `equality`, and `inequality` to the function `dense_sdp_relax` to generate SDP relaxations. We recommend the user to run the script [`example_bqp.m`](https://github.com/MIT-SPARK/CertifiablyRobustPerception/blob/master/BinaryQuadraticProgram/example_bqp.m) to check how to solve the `SDP` data using MOSEK, and to see that SDP relaxations can actually solve BQP problems to global optimality.

### Sparse Relaxation (Basis Reduction)
Coming soon.
