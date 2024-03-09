# DR-Regret
====================================================

This repository includes the MATLAB codes to perform the simulations in **[Wasserstein Distributionally Robust Regret Minimization][paper_link]**.

## 1. Requirements
To run our codes, the following softwares must be installed:
- **[MATLAB][MATLAB]**
- **[CPLEX Optimization Studio][CPLEX]**

which are free for students and academics. Regarding their versions, we used R2023a and 12.10, respectively. 

## 2. Simulations 
The folder has seven script files. By running `main.m`, which contains all the other six, you can obtain solutions to the three decisioin-making models and compute, for each solution, the true ex-ante regret as well as the upper and lower bounds on the worst-case ex-ante regret.
Specifically, `solution_dro_exantereg.m`, `solution_dro_expostreg.m` and `solution_dro_cost.m` return solutions to the proposed model, M1 and M2, respectively.
Executing `compute_maxreg_ub.m` and `compute_maxreg_lb.m` yield the upper and lower bounds on the worst-case ex-ante regret. 
The problem parameters such as the mean of uncertain demand can be changed by modifying `problem_newsvendor.m`.

[paper_link]: ..
[MATLAB]: https://matlab.mathworks.com
[CPLEX]: https://www.ibm.com/products/ilog-cplex-optimization-studio
