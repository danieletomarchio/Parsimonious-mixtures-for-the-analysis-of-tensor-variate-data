# Parsimonious-mixtures-for-the-analysis-of-tensor-variate-data

This repository contains the code for fitting parsimonious mixtures of tensor-variate shifted exponential normal (TVSEN), tensor-variate tail inflated normal (TVTIN), and tensor-variate normal (TVN) distributions.
In the following, you find a description of the functions (and their arguments) contained in the Tensor-mixtures.R file.

## TensMixtPars_init ##

### Description ###

Runs the initialization of the EM-based algorithms used for fitting parsimonious mixtures of TVSEN, TVTIN, and TVN distributions. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

TensMixtPars_init (X, k, mod.list = "all", nstartR = 50, nThreads = 1, density = "TVN")

### Arguments ###

* X: An array with at least D=4 dimensions, where the last is occupied by the statistical units. For example, a 4-way array has D1 x D2 x D3 x D4 dimensions, where D1, D2, and D3 are reserved for the variables, whereas D4 considers the statistical units.
* k: A vector (or a number) containing the groups to be tried.
* mod.list: A list, containing the specific parsimonious models to be initialized, or the character "all" used to consider all the possible parsimonious combinations.
When the list is used, each of its elements must be a character vector indicating the specific parsimonious structure. Possible values for each of the first D-1 elements of the character vectors are: "II", "EI", "VI", "EE", "VE", "EV", "VV". The Dth element of the character vectors can have the following possible values: "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV".
* nstartR: An integer specifying the number of random starts to be considered. Default value is 50. 
* nThreads: A positive integer indicating the number of cores used for running in parallel.  
* density: A character specifying the tensor-variate distribution to consider. Possible values are: "TVN", "TVSEN", and "TVTIN".

## TensMixtPars_fit ##

### Description ###

Fits, by using EM-based algorithms, parsimonious mixtures of TVSEN, TVTIN, and TVN distributions. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

TensMixtPars_fit (X, init.par = NULL, tol = 0.001, maxit = 500, nThreads = 1)

### Arguments ###

* X: An array with at least D=4 dimensions, where the last is occupied by the statistical units. For example, a 4-way array has D1 x D2 x D3 x D4 dimensions, where D1, D2, and D3 are reserved for the variables, whereas D4 considers the statistical units.
* init.par: The output of the TensMixtPars_init() function.
* tol: Threshold for EM-based algorithms convergence. Default value is 0.001.
* maxit: Maximum number of iterations for the EM-based algorithms. Default value is 500.
* nThreads: A positive integer indicating the number of cores used for running in parallel.

## extract.bestM ##

### Description ###

This function extracts the top/best fitting models according to the Bayesian information criterion (BIC).

### Usage ###

extract.bestM (results, top = 1, short = FALSE)

### Arguments ###

* results: The output of the TensMixtPars_fit() function.
* top: A number indicating how many models to extract from the ranking provided by the BIC. 
* short: A logical indicating whether the extraction refers to the models fitted until convergence (FALSE) or whether the extraction refers to the models fitted up to a short number of iterations by using our fitting strategy (TRUE). Default value is FALSE.

## extract.shortEM ##

### Description ###

This function extracts the names and k of the best fitting models provided by the extract.bestM function.

### Usage ###

extract.shortEM (Extract.bestM)

### Arguments ###

* Extract.bestM: The output of the extract.shortEM() function.

## filt.init ##

### Description ###

This function extracts only the initializations to be used in the second step of our fitting strategy, according to the results provided by the extract.shortEM() function.

### Usage ###

filt.init (res.init, Extract.shortEM)

### Arguments ###

* res.init: The output of the TensMixtPars_init() function.
* Extract.shortEM: The output of the extract.shortEM() function.

## Trnorm ##

### Description ###

This function generates random observations from a TVN distribution.

### Usage ###

Trnorm (n, m.vec, mu = array(0, m.vec), Sigma.list = NULL)

### Arguments ###

* n: The number of statistical units.
* m.vec: A vector of at least three elements specifying the dimension of each variable.
* mu: An array of at least three elements specifying the mean matrix of each dimension.
* Sigma.list: A list of at least three elements specifying the covariance matrix of each dimension.

## TrSenTin ##

### Description ###

This function generates random observations from a TVSEN and TVTIN distributions.

### Usage ###

TrSenTin (n, M, Sigma, theta, density)

### Arguments ###

* n: The number of statistical units.
* M: An array of at least three elements specifying the mean matrix of each dimension.
* theta: A number referring to the tailedness parameter of the chosen distribution.
* density: A character specifying the tensor-variate distribution to consider for generating the data. Possible values are: "TVSEN" and "TVTIN". 

# Example files

The example files allow to simulate and fit our models.
