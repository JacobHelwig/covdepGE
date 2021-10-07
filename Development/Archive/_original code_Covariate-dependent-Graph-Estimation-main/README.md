# Covariate-dependent-Graph-Estimation

There are two demo codes in this repository which showcases the practical performance of the graph estimator proposed in our paper "An approximate Bayesian approach to covariate dependent graphical modeling". The code discrete_covariate_demo.R considers the case of discrete covariates. The code cont_covariate_attempt2.R considers the toy example presented in the paper with continuous covariates. The cov_vsvb.R and ELBO_calculator are functions called by the demo codes. Specifically, the function cov_vsvb updates the variational parameters and returns the converged estimates for a single graph.

Overview of discrete_covariate_demo.R:
=====================================
One can simply run the demo as is to get some demo examples and some visual results through a heatmap and histograms.

## Data generation
In this file, it is assumed that there are 2 discrete covariate levels. The data are generated from two different covariance matrices as an example, controlled by a ![equation](https://latex.codecogs.com/gif.latex?%5Clambda) parameter. Depending on whether ![equation](https://latex.codecogs.com/gif.latex?%5Clambda_1%3D%5Clambda_2), we have the covariate independent model or the covariate dependent model. Set no. of subjects in study to be `n` and number of variables to be `p+1`.

#1. Covariate independent model

![equation](https://latex.codecogs.com/gif.latex?%5Clambda_1%3D%5Clambda_2%3D%7B%5B15%7B%5Cbf%201%7D_4%7E%7E%20%2C%20%7E%7E%7B%5Cbf%200%7D_%7Bp-3%7D%5D%7D%5E%7B%5Cmathrm%7BT%7D%7D)

#2. Covariate dependent model

![equation](https://latex.codecogs.com/gif.latex?%5Clambda_%7B1%7D%3D%20%7B%5B15%7B%5Cbf%201%7D_4%7E%7E%20%2C%20%7E%7E%7B%5Cbf%200%7D_%7Bp-3%7D%5D%7D%5E%7B%5Cmathrm%7BT%7D%7D)

![equation](https://latex.codecogs.com/gif.latex?%5Clambda_%7B2%7D%3D%20%7B%5B%7B%5Cbf%200%7D_%7Bp-3%7D%7E%7E%20%2C%20%7E%7E15%7B%5Cbf%201%7D_%7B4%7D%5D%7D%5E%7B%5Cmathrm%7BT%7D%7D)

The precision matrix ![equation](https://latex.codecogs.com/gif.latex?%5COmega_1%20%3D%20%5Clambda_1%5Clambda_1%5ET%20&plus;%2010%20%5Cmathbb%7BI%7D). Similarly, ![equation](https://latex.codecogs.com/gif.latex?%5COmega_2%20%3D%20%5Clambda_2%5Clambda_2%5ET%20&plus;%2010%20%5Cmathbb%7BI%7D).
Let ![equation](https://latex.codecogs.com/gif.latex?%5CSigma_1%3D%20%5COmega_1%20%5E%7B-1%7D), and ![equation](https://latex.codecogs.com/gif.latex?%5CSigma_2%3D%20%5COmega_2%20%5E%7B-1%7D).
We generate `n/2` samples from ![equation](https://latex.codecogs.com/gif.latex?%5Cmathcal%7BN%7D%280%2C%5CSigma_1%29) and `n/2` samples from ![equation](https://latex.codecogs.com/gif.latex?%5Cmathcal%7BN%7D%280%2C%5CSigma_2%29) to form our dataset.

## Generating the covariates
We generate an ![equation](https://latex.codecogs.com/gif.latex?n%20%5Ctimes%20p) covariate matrix with entries = `-0.1` for each of the p variables for the population with precision matrix ![equation](https://latex.codecogs.com/gif.latex?%5COmega_1) and entries = `0.1` for each of the p variables for the population with precision matrix ![equation](https://latex.codecogs.com/gif.latex?%5COmega_2). Thus for each variable, the covariate value for an individual is univariate. As an example, the covariate attached to the FOXC2 protein expression of patient 1 is the univariate FOXC2 RNA expression for the same patient. If instead we used both the RNA expression and CNV expression for FOXC2 gene for the same patient, we would have a two-dimensional covariate attached to the data.

## Overview of the algorithm
1. Fix a variable  `j` as response, and the remaining `p` variables as predictor. (Recall there are `p+1` variables total.
2. From the covariate matrix, define an ![equation](https://latex.codecogs.com/gif.latex?n%20%5Ctimes%20n) weight matrix where the `i`th row describes the weight vector associated with the n subjects relative to subject `i`. The weights for this model are chosen with an ad-hoc bandwidth value of 0.1. Technically, one can perform a density estimation on the covariate space, but since its basically discrete, we choose a small value of 0.1 
3. Choose the hyperparameter values ![equation](https://latex.codecogs.com/gif.latex?%5Cpi) and ![equation](https://latex.codecogs.com/gif.latex?%5Csigma%5E2) following Carbonetto Stephens and ![equation](https://latex.codecogs.com/gif.latex?%5Csigma_%7B%5Cbeta%7D%5E2) over a grid.
4. Call the cov_vsvb function to update the following variational parameters : The ![equation](https://latex.codecogs.com/gif.latex?n%20%5Ctimes%20p) matrices alpha, mu and S_sq where the `i`th row corresponds to the inclusion probability of the `p-1` predictor variables, mean and standard deviation for the `i`th subject.  
5. Loop over the `p` variables as response to get the ![equation](https://latex.codecogs.com/gif.latex?%28p&plus;1%29%20%5Ctimes%20p) matrices corresponding to each of the `n` subjects in the study.
6. Assume that the diagonal elements in the inclusion probability matrices for each individual is `0`, and apply the post processing ![equation](https://latex.codecogs.com/gif.latex?%5Calpha_%7Bjk%7D%5E*%3D%28%5Calpha_%7Bjk%7D%20&plus;%20%5Calpha_%7Bkj%7D%29/2) to symmetrize the matrix.
7. Set the dependence graph to be ![equation](https://latex.codecogs.com/gif.latex?%5Cmathbb%7BI%7D%5C%7B%5Calpha_%7Bjk%7D%5E*%3E0.5%5C%7D).


## Remarks 
The code calls the cov_vsvb function, which updates the variational parameters and returns the final values of the estimates corresponding to a single graph which corresponds to a fixed individual in the study. By going through a loop, one can calculate the graph estimates corresponding to every individual in the study, and can also be parallelized since the updates are independent. However, the parallelization is yet to be implemented.
The cov_vsvb function itself calls the ELBO_calculator function, which calculates the ELBO corresponding to the current values of the variational parameters for a specific graph corresponding to a single individual. Note the contribution of every individual in the study to the ELBO of the parameters corresponding to a single individual, facilitating the borrowing of information.


Overview of cont_covariate_demo.R:
====================================
In this file, instead of discrete covariate values, we have three clusters of covariate values, and the individuals in the study have covariate values belonging to one of the three clusters. 

## Data generation
We set `n=180` and `p=4`, i.e. there are `5` variables in this toy example. We have univariate covariates associated with every individual. 

`Z ~ Uniform (-1,-0.3) U (-0.23, 0.33) U (0.43, 1)` corresponding to three well-separated clusters. 
We have the precision matrix defined as a function of the covariate Z, through the var_cont function as 

![equation](https://latex.codecogs.com/gif.latex?%5COmega%5Ei%28z_i%29%3D%20%5Cbegin%7Bbmatrix%7D%202%20%26%20%5Comega_%7B1%2C2%7D%28z_i%29%20%26%20%5Comega_%7B1%2C3%7D%28z_i%29%20%260%20%260%20%5C%5C%20%5Comega_%7B1%2C2%7D%28z_i%29%20%26%202%20%26%201%20%260%20%260%5C%5C%20%5Comega_%7B1%2C3%7D%28z_i%29%20%26%201%20%26%202%20%26%200%20%26%200%5C%5C%200%20%26%200%20%260%20%262%20%260%5C%5C%200%26%200%26%200%20%26%200%20%26%202%20%5Cend%7Bbmatrix%7D)

where

![equation](https://latex.codecogs.com/gif.latex?%5Comega_%7B1%2C2%7D%28z_i%29%3D%5Cbegin%7Bcases%7D%20%241%24%20%26%20%5Ctext%7B%20if%20%24-1%20%3C%20z_i%20%3C-0.33%24%7D%20%5C%5C%201-%20%5Cfrac%7Bz_i%20&plus;%200.23%7D%7B0.33&plus;0.23%7D%20%26%20%5Ctext%7B%20if%20%24-0.23%3Cz_i%3C0.33%24%7D%5C%5C%20%240%24%20%26%20%5Ctext%7B%20if%20%240.33%3Cz%3C1%24%7D%20%5Cend%7Bcases%7D%20%2C%20%5Comega_%7B1%2C3%7D%28z_i%29%3D%5Cbegin%7Bcases%7D%20%240%24%20%26%20%5Ctext%7B%20if%20%24-1%20%3C%20z_i%20%3C-0.33%24%7D%20%5C%5C%20%5Cfrac%7Bz_i%20&plus;%200.23%7D%7B0.33&plus;0.23%7D%20%26%20%5Ctext%7B%20if%20%24-0.23%3Cz_i%3C0.33%24%7D%5C%5C%20%241%24%20%26%20%5Ctext%7B%20if%20%240.33%3Cz%3C1%24%7D.%20%5Cend%7Bcases%7D)

The rest of the steps are identical to the discrete covariate case.
