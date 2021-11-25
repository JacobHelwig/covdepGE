# `covdepGE:` Covariate Dependent Graph Estimation

## Installation

Run the following in `R`:

`devtools::install_github("JacobHelwig/covdepGE")`

## Overview

Suppose **X** ∈ ℝ<sup>*n* × *p*</sup> is a data matrix of independent
observations **X** = (**x**<sub>1</sub>,...,**x**<sub>*p*</sub>), where,
for *j* ∈ 1, ..., *p*:

**x**<sub>*j*</sub> ∼ 𝒩(*μ*<sub>*j*</sub>,*Σ*<sub>*j*, *j*</sub>), **x**<sub>*j*</sub> ∈ ℝ<sup>*n*</sup>      **X** ∼ 𝒩(*μ*,*Σ*)

The conditional dependence structure of
**x**<sub>**1**</sub>, ..., **x**<sub>**p**</sub> can be modeled as an
undirected graph 𝒢 such that:

$$\\mathcal G\_{i,j} = 
\\begin{cases}
    1 & \\iff \\mathrm{Cov}(\\pmb{x_i}, \\pmb{x_j})\\neq 0 
    \\\\
    0 & \\textrm{otherwise} 
\\end{cases}$$

That is, there is an edge between the **x**<sub>**i**</sub> and
**x**<sub>**j**</sub> nodes if, and only if, these variables are
dependent on each other given all other variables.

Let **Z** be an *n* × *p*′ matrix of extraneous covariates. Further
suppose that the conditional dependence structure of **X** is not
homogeneous across the individuals, and is instead a continuous function
of the extraneous covariates **Z**(1). Then, this methodology aims to
estimate a graph for each of the individuals, possibly unique to the
individual, such that similar estimates are made for those who are
similar to one another in terms of the extraneous covariates.

For an example application, see (1), wherein the sample was composed of
healthy and cancerous individuals,
**x**<sub>1</sub>, ..., **x**<sub>8</sub> were protein expression levels
of 8 genes, and **Z** was the copy number variation of a gene **z**
associated with cancer,
**z** ∉ {**x**<sub>1</sub>, ..., **x**<sub>8</sub>}.

## Functionality

The main function, `covdepGE::covdepGE(`**X**, **Z**`)`, estimates the
posterior distribution of the graphical structure 𝒢<sub>*l*</sub> for
each of the *n* individuals using a variational mean-field
approximation. The function will output *n* *p* × *p* symmetric matrices
𝒜<sub>*l*</sub>, where 𝒜<sub>*i*, *j*</sub><sup>(*l*)</sup> is the
posterior inclusion probability of an edge between the node representing
the *i*-th variable and the node representing the *j*-th variable.

## To-do

-   Complete documentation

-   Add compatibility checks

-   Add warnings for when the optimal `sigmabeta_sq` OR `pi_est` is
    chosen at either of the grid end points; a `warning` argument for
    when users would like to toggle off the warnings. Should only be
    enabled when either `sigmabeta_sq` or `pi_est` is longer than 2.

-   Create a vignette demonstrating usage on a simple simulated dataset.

-   Change the `idmod` probs to logbase 10

-   Remove `CS` argument

-   Parallelization of the “main loop” over the predictors in
    `covdepGE_main.R`). This is complicated by the `C++` code, however,
    two potential solutions are:

    -   <span style="color: blu">[StackOverflow
        suggestion](https://stackoverflow.com/questions/69789634/parallelization-of-rcpp-without-inline-creating-a-local-package?noredirect=1#comment123649680_69789634)</span>

    -   <span
        style="color: blu">[RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)</span>

    This is a finishing touch and **likely will not be implemented in
    the package by the end of the semester.**

<div class="thebibliography">

1 (1) Dasgupta, Sutanoy, et al. “An approximate Bayesian approach to
covariate dependent graphical modeling." 2021

</div>
