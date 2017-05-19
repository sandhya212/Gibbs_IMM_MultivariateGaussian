Readme
=======
Source code that implements an Infinite Mixture Model for multivariate Gaussian data and infers the mixture model parameters via Gibbs MCMC. 

Software provided in R.

**Dependencies**

1. R - install R from https://cran.r-project.org/. A version higher than 2.12.0 is recommended.
sudo yum install R
2. Once installed, open a terminal and at the command prompt, type R. 
3. At the R prompt: Install the following R packages by issuing command:

>install.packages(c("MCMCpack","mvtnorm","coda","Matrix","Rtsne","lattice","MASS","bayesm","robustbase","chron","mnormt","schoolmath","devtools"))


> **Note:**
Until here, this is a one-time activity. 


**How to run code**


1. Clone/Download code repository

2. At the R prompt: issue command
>rm(list=ls())
graphics.off()

3. Issue setwd() to point to the path where the code repository resides. 
eg. if the code is downloaded at “/User/Downloads/IMM_Code/“, then type
>working_path <- "/Users/Downloads/IMM_Code/“; 
setwd(working_path);

4. At the R prompt: Issue command
> source(start_file.R)

**Output**
______

1. ~/output/plots/ folder gives the various plots
3. In the R terminal, the latent variables of interest can be obtained by issuing:

>- Class               ##  *(class assignment of the cells)*
>- mu                         ##  *(inferred means. dim(mu_final) <- numgenes x K)*
>- Sigma                    ##  *(inferred covariances)*
