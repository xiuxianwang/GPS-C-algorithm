The methods used to estimate the Gaussian process parameters, in the SKfit.m files, are different for homoscedastic and heteroscedastic problems. 

For homoscedastic problems, the MLE method used in the SKO algorithm of Huang et al. (2006) is used to estimate the Gaussian process parameters (mu_0, tau^2 and theta) and the unknown simualtion noise variance simultaneously.

For heteroscedastic problems, a kernel-based method is used first to estimate the variances of simulation noises (provided in Section 5.2.1). Then, based on these estimated variances, the Gaussian process parameters are estimated using a different MLE method. This MLE method is used in the stochastic kriging of Ankenman et al. (2010).

The methods used to find the current best solution of each iteration are also different depending on the dimensions of test problems. 

In two-dimensional problems, a dense grid is constructed in the feasible regions. In each iteration, the grid-best solution is identified first, then the Matlab solve fmincon is used to further improve this grid-best solution in its neighborhood.

In high-dimensional problems, the Matlab solve fmincon is used directly to search good solutions in the whole feasible region. 