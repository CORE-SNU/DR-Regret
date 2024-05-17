# DR-Regret
====================================================

This repository includes the MATLAB codes to perform the simulations in **[Wasserstein Distributionally Robust Regret Minimization][paper_link]**.

## 1. Requirements
To run our codes, the following softwares must be installed:
- **[MATLAB][MATLAB]**
- **[CPLEX Optimization Studio][CPLEX]**

which are free for students and academics. Regarding their versions, we used R2023a and 12.10, respectively. 

## 2. Simulations 
The folder has seven script files. By running main.m, which contains all the other six, you can obtain solutions to the three decisioin-making models and compute the true ex-ante regret as well as the upper and lower bounds on the worst-case ex-ante regret. Specifically, main.m and problem_newvendor.m are described as follows: 
- main.m:
  - `line 3-9`: parameters related to the newsvendor problem and the algorithm such as the sample size, the set of Wassersetin ball's radii, and the convergence tolerance are set.
  - `line 12-15`: sample indices are randomly selected.
  - `line 19-31`: solution_dro_exantereg.m is run, by which solutions and average upper bounds on the maximum ex-ante regret of our model are computed.
  - `line 34-43`: solution_dro_expostreg.m is run, by which solutions of M2 are computed.
  - `line 46-61`: solution_dro_cost.m is run, by which solutions of M1 are computed.
  - `line 65-76`: compute_maxreg_ub.m is run, by which average upper bounds on the maximum ex-ante regret of M2 are computed.
  - `line 79-90`: compute_maxreg_ub.m is run, by which average upper bounds on the maximum ex-ante regret of M1 are computed.
  - `line 94-105`: compute_maxreg_lb.m` is run, by which average lower bounds on the maximum ex-ante regret of our model are computed.
  - `line 108-119`: compute_maxreg_lb.m` is run, by which average lower bounds on the maximum ex-ante regret of M2 are computed.
  - `line 122-133`: compute_maxreg_lb.m` is run, by which average lower bounds on the maximum ex-ante regret of M1 are computed.
  - `line 143-161`: average true ex-ante regret of our model are computed.
  - `line 164-182`: average true ex-ante regret of M2 are computed.
  - `line 185-203`: average true ex-ante regret of M1 are computed.
- problem_newsvendor.m:
  - `line 23-35`: the cost function $f(x,\xi)$ is represented as the pointwise maximum of affine functions to compute solutions to M1.
  - `line 38-43`: the true distribution of uncertainty is modeled.
  - `line 46-48`: a discrete distribution approximating the true distribution of uncertainty is obtained.
  - `line 51-93`: the minimum expected cost of the newsvendor problem for the true distribution of uncertainty is computed using the alternating direction method of multipliers (ADMM).


[paper_link]: ..
[MATLAB]: https://matlab.mathworks.com
[CPLEX]: https://www.ibm.com/products/ilog-cplex-optimization-studio
