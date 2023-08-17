In directories "rate_two dimensional" and "rate_three dimensional", the test problems are generated from Gaussian processes. Therefore, some additional or different files are used.

Files SKfit_create_sample_path.m and SKpredict_sample_path.m are used to generate sample paths by calculating the conditional mean values based on observations collected from a multivariate normal distribution.

File fun_loc_many.m is used to generate noisy observations from the generated sample paths.