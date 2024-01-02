# sampling large hyperplane truncated MVN distrbution using Karhunen-Loève expansion (KLE) and Matheron's update rule (MUR)
This repository includes R code for sampling from a high-dimensional hyperplane-truncated multivariate normal distribution using the Karhunen-Loève expansion (KLE) and Matheron's update rule (MUR). Additionally, it provides the "LS.KLE" function for generating a very large Gaussian vector extracted from a stationary Gaussian process.

# References
Maatouk, H. and Rullière, D. and Bay, X. (2023). "Sampling large hyperplane‐truncated multivariate normal distributions". Computational Statistics. [doi](https://link.springer.com/article/10.1007/s00180-023-01416-7)

Wood, A. and Chan, G. (1994). "Simulation of Stationary Gaussian Processes in [0,1]^d". Journal of Computational and Graphical Statistics. [doi](https://www.jstor.org/stable/1390903)

# Images
[MyImage](https://github.com/maatouk/sampling_large_hyperplane_tMVN/blob/main/Matern3split5sim-eps-converted-to.pdf): GP sample paths when the domain is split into three subdomains. Dashed (resp. solid) curves represent the five paths after (resp. before) conditioning. 

[MyImage](https://github.com/maatouk/sampling_large_hyperplane_tMVN/blob/main/FFTvsLSKLE-eps-converted-to.pdf): Average time of sampling a MVN over 25 replicates as a function of the dimension N. The proposed approach, LS.KLE, has been compared to the Fast Fourier Transform (FFT) developed in Wood and Chan (1994) 


# Description of the associated R files:
1. 'all_models.R' file contains two models that implement the large hyperplane-truncated multivariate normal distribution with and without hyperparameters updates.
2. 'Illustr_ex_COSTA.R' file includes all synthetic and real data examples featured in Maatouk et al. (2023).
3. 'all-fcts.R' file encompasses all the fundamental functions required to implement our approach. It includes:
   
i. The naive Matheron's update rule (MUR).

ii. The naive Karhunen-Loève Expansion (KLE).

iii. The proposed LS.KLE for drawing very large Gaussian priors.

iv. nu.MH2 for drawing posterior samples of the hyperparameters using the Metropolis–Hastings method and inverse Cholesky factor.

v. LS.KLE_MUR for sampling from large hyperplane-truncated multivariate normal distributions.

vi. samp.WC for generating Gaussian priors using the Fast Fourier Transform method of Wood and Chan (1994).

   For more details on the codes or the functions, refer to the associated R files.
