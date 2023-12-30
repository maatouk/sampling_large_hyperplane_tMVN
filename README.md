# sampling large hyperplane truncated MVN distrbution using Karhunen-Loève expansion (KLE) and Matheron's update rule (MUR)
This repository includes R code for sampling from a high-dimensional hyperplane-truncated multivariate normal distribution using the Karhunen-Loève expansion (KLE) and Matheron's update rule (MUR). Additionally, it provides the "LS.KLE" function for generating a very large Gaussian vector extracted from a stationary Gaussian process.

# References
Maatouk, H. and Rullière, D. and Bay, X. (2023). "Sampling large hyperplane‐truncated multivariate normal distributions", To appear in Computational Statistics. [doi](https://link.springer.com/article/10.1007/s00180-023-01416-7)

Wood, A. and Chan, G. (1994). "Simulation of Stationary Gaussian Processes in [0,1]^d". Journal of Computational and Graphical Statistics. [doi](https://www.jstor.org/stable/1390903)


# Description of the associated R files:
1. 'all_models.R' file contains two models that implement the large hyperplane-truncated multivariate normal distribution with and without hyperparameters updates.
2. 'Illustr_ex_COSTA.R' file contains all the numerical examples provided in Maatouk et al. (2023).
3. 'all-fcts.R' file contains all the base functions need to implemented our approach, including the naive Matheron's update rul (MUR), the naive KLE, the proposed LS.KLE (draws very large Gaussian prior), nu.MH2 (draws posterior samples on the hyperparameters using Metropolis–Hastings and inverse Cholesky factor), LS.KLE_MUR (draws large hyperplane-truncated multivariate normal distribution), samp.WC (draws Gaussian prior using the Fast Fourier Transform of Wood and Chan [1994]).

   For more details on the codes or the functions, refer to the associated R files.
