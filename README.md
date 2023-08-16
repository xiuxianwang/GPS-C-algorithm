# GPS-C-algorithm

On this site, we include all the codes used to conduct the numerical experiments of the paper entitled "Gaussian Process based Random Search for Continuous Optimization vis Simulation". 
In Section 6 of the paper, we conduct a series of numerical experiments to show the empircal rate of convergence and finite-sample performance of the GPS-C algorithm. All codes are 
written in Matlab. 

Directory "GPS-C" contains the codes of the GPS-C algorithm that are used in Sections 6.2.1 and 6.2.2 to test its empirical performance under both homoscedastic and heteroscedastic contexts.

Directory "Rate" contains the codes of GPS-C algorithm that are used in Sections 6.1 and 6.3 to test its empirical rate of convergence for solving low-dimensional (2 and 3) generated 
problems (randomly generated from Gaussian process) and the Weighted-Sphere problems of different (but higher) dimensions.

Directory "Other algorithms" contains all the codes of the SKO algorithm, the KGCP algorithm, the ASR algorithm and the SOSA algorithm that are used in Section 6.2 to provide benchmarks 
for the GPS-C algorithm.


For the codes of the GPS-C algorithm, CGPS.m is the main file. Files Fun_xxx.m are the functions of test problems with simulation noises. Files Fun_xxx_free are the functions of test problems 
without simulation noises, which are used to evaluate the quality of obtained solutions.

Files SKfit.m and SKfit2.m are used to estimate the Gaussian process parameters required by the GPS-C algorithm. File Skpredict.m is used to predict the conditional mean values of each sampled 
solution. File CalculateMSE.m is used to calculate the conditional variance of each sampled solution. This file is used in the sampling component. All these codes are revised based on the codes of 
stochstic kriging of Ankenman et. al. (2010). Files LogPL.m, corrcubR.m and correxpR.m are the functions required by Files SKfit.m and SKpredict.m.

Files Matrix_inverse.m and Matrix_inverse_1.m are used to calculate and record the inverse of the covariance matrix of each iteration. File Matrix_inverse.m use the inverse function of Matlab _inv()_ 
directly. File Matrix_inverse_1.m uses block matrix inversion to reduce computational complexity.

Files SKpredict_fmincon.m and find_current_best.m are used to find the current optimal solution of the contructed Gaussian-surrogate model. File SKpredict_fmincon.m is used in two-dimensional problems, where 
it is possible to construct and evaluate a dense grid, to find current best solutions around the grid-best solution. File find_current_best.m is used in high-dimensional problems to search good solutions 
in the whole feasible region. In both files, the Matlab optimizer _fmincon()_ is used.

File Sampling.m is used to sample solutions in each iteration. The Markov Chain Coordinate Sampling (MCCS) scheme is adopted in this paper.
