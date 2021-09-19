# Certifiable Outlier-Robust Geometric Perception

## About
This repository holds the implementation for certifiably solving outlier-robust geometric perception problems to global optimality. The certifiable outlier-robust geometric perception framework contains two main modules:

- A sparse semidefinite programming relaxation (SSR) scheme that relaxes nonconvex outlier-robust perception problems into convex semidefinite programs (SDPs); and

- A novel SDP solver, called [STRIDE](https://github.com/MIT-SPARK/STRIDE), that solves the generated SDPs at an unprecedented scale and accuracy.

We proposed a preliminary version of the SSR scheme in our [NeurIPS 2020 paper](https://arxiv.org/abs/2006.06769), and released a certifier (that certifies if a given estimate is optimal) based on Douglas-Rachford Splitting (DRS). Please switch to the `NeurIPS2020` branch of this repo to checkout the NeurIPS2020 implementation.

If you find this library helpful or use it in your projects, please cite:
```bibtex
@article{Yang21arXiv-certifiableperception,
  title={Certifiable Outlier-Robust Geometric Perception: Exact Semidefinite Relaxations and Scalable Global Optimization},
  author={Yang, Heng and Carlone, Luca},
  journal={arXiv preprint arXiv:2109.03349},
  year={2021}
}
```

## Tutorial on Semidefinite Relaxation
A general polynomial optimization problem (POP) is an optimization problem of the following standard form
```
minimize_{x \in R^d}    f(x)
subject to              h_i(x) =  0, i=1,...,l_h
                        g_j(x) >= 0, j=1,...,l_g
```
where `x \in R^d` is a `d`-dimensional decision variable, `f(x)` is a scalar polynomial objective function, `h_i(x),i=1,...,l_h` are scalar polynomial functions that define equality constraints, and `g_j(x),j=1,...,l_g` are polynomial functions that define inequality constraints. POPs are in general NP-hard problems. For example, one can easily see that binary quadratic programming is an instance of POP, because binary constraints `x_i \in {1, -1}, i=1,...,d` can be easily written as polynomial equalities `h_i(x) = x_i^2 - 1 = 0,i=1,...,d`. 

### Dense Relaxation (Lasserre's Hierarchy)

### Sparse Relaxation (Basis Reduction)
