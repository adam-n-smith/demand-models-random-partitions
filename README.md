# Demand Models with Random Partitions

This repository contains the R and Rcpp source code to implement the location-scale partitioning (LSP) methodology of Smith and Allenby (2019).

## List of Files 

- `dmrpfunctions.cpp`

  This is an Rcpp file containing four sections of code. A description of each section along with its exported functions is provided below.

  **Section 1** contains core functions used in other parts of the file. 
  - `compsm`: compute pairwise similarity matrix given draws from a partition distribution

  **Section 2** contains functions to sample from and evaluate PMFs for the DP, LSP, LSPx, and block LSP models.
  - `rdp`: sample from the DP partition distribution
  - `rlsp`: sample from the LSP distribution
  - `rlspx`: sample from the LSPx distribution
  - `rblocklsp`: sample from the LSP distribution at given starting/ending point
  - `ddp`: evaluate PMF of the DP partition distribution
  - `dlsp`: evaluate PMF of the LSP distribution
  - `dlspx`: evaluate PMF of the LSPx distribution
  - `dblocklspx`: evaluate PMF of the LSP distribution at given starting/ending point

  **Section 3** contains functions to carry out a simulation study to compare LSP proposals with alternative proposal mechanisms.
  - `regll`: evaluates the log likelihood of a univariate regression model with a nonlinear mean function
  - `drawblocks`: draws block structure used for block LSP proposals
  - `reglsp`: MCMC sampler for regression model using LSP proposals
  - `regblocklsp`: MCMC sampler for regression model using block LSP proposals

  **Section 4** contains the MCMC samplers for an isolated demand model with a random partition and various partitioning priors.
  - `mvregll`: evaluates the log likelihood of a multivariate regression model
  - `mvregmcmc`: MCMC sampler for an unrestricted multivariate regression model
  - `isomvregdpmcmc`: MCMC sampler for an isolated demand model with a random partition and DP prior
  - `isomvreglspmcmc`: MCMC sampler for an isolated demand model with a random partition and LSP prior
  - `isomvreglspxmcmc`: MCMC sampler for an isolated demand model with a random partition and LSPx prior

- `LSPproperties.R`

  This is an R file that illustrates the properties of the LSP and LSPx models (figures 1 and 2). There is also code to compare the LSP distribution to the DP, ddCRP, and EPA partition distributions (figure 3).

- `compareproposals.R`
 
  This is an R file that compares the performance of MCMC algorithms using LSP and block LSP proposals with split-merge proposals and Gibbs proposals. Simulation studies are conducted for both small n (figure 6) and large n (figure 7).

- `simulationstudy.R`

  This is an R file that creates a simulated data set and then fits an unrestricted demand model and an isolated demand model (with both DP and LSP priors) to the data. There is code to check the recovery of all model parameters, plot posterior pairwise similarity matrices, and compute model fit statistics.

  Note: To get an idea of computation time, if a data set is generated with n=20, K=5, and T=100, then the isomvreglspmcmc sampler generates roughly 1,000 posterior draws per second. Computation time will increase as n increases, K decreases, or T increases. 


## References
- Adam N. Smith & Greg M. Allenby (2019) Demand Models With Random Partitions, *Journal of the American Statistical Association*, [DOI: 10.1080/01621459.2019.1604360](https://doi.org/10.1080/01621459.2019.1604360)
