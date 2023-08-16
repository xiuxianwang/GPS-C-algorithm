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

Files SKfit.m and SKfit_2.m are used to estimate the Gaussian process parameters required by the GPS-C algorithm. SKfit_2.m is used in 
