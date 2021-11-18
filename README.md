# `covdepGE:` Covariate Dependent Graph Estimation

## Overview

Suppose $\pmb{X}\in \mathbb{R}^{n\times p}$ is a data matrix of
independent observations
$\pmb{X}= \left(\pmb{x}_1,..., \pmb{x}_p\right)$, where, for
$j\in1,...,p$:

$$\pmb{x}_j\sim\mathcal{N}\left(\mu_j,\Sigma_{j,j}\right), \pmb x_j\in \mathbb{R}^n\hspace{2cm}\pmb{X}\sim\mathcal{N}\left(\mu,\Sigma\right)$$

Let $\pmb{Z}$ be an $n\times p'$ matrix of extraneous covariates. The
conditional dependence structure of $\pmb{x_1},...,\pmb{x_p}$ can be
modeled as an undirected graph $\mathcal{G}$ such that:

$$\mathcal G_{i,j} = 
\begin{cases}
    1 & \iff \mathrm{Cov}(\pmb{x_i}, \pmb{x_j})\neq 0 
    \\
    0 & \textrm{otherwise} 
\end{cases}$$

That is, there is an edge between the $\pmb{x_i}$ and $\pmb{x_j}$ nodes
if, and only if, these variables are dependent on each other given all
other variables.

Further suppose that the conditional dependence structure of $\pmb{X}$
is not homogeneous across the individuals, and is instead a continuous
function of the extraneous covariates $\pmb Z$[@covDepGM]. Then, this
methodology aims to estimate a graph for each of the individuals,
possibly unique to the individual, such that similar estimates are made
for those who are similar to one another in terms of the extraneous
covariates.

For an example application, see [@covDepGM], wherein the sample was
composed of healthy and cancerous individuals, $\pmb x_1,...,\pmb x_8$
were protein expression levels of 8 genes, and $\pmb Z$ was the copy
number variation of a gene $\pmb z$ associated with cancer,
$\pmb z\not\in \{\pmb x_1,...,\pmb x_8\}$.

## Functionality

The main function, `covdepGE::covdepGE(`$\pmb X,\pmb Z$`)`, estimates
the posterior distribution of the graphical structure $\mathcal G_l$ for
each of the $n$ individuals using a variational mean-field
approximation. The function will output $n$ $p\times p$ symmetric
matrices $\mathcal{A}_l$, where ${\mathcal{A}_{i,j}^{(l)}}$ is the
posterior inclusion probability of an edge between the node representing
the $i$-th variable and the node representing the $j$-th variable.

## To-do

-   Implement KDE to obtain individual specific bandwidths. Currently,
    only one global bandwidth is used, and it requires that the user
    specify it (or use the default, deterministic value). This
    modification will estimate a unique bandwidth for each individual
    dependent on the empirical density of $\pmb Z$ resulting from KDE.

-   Create a vignette demonstrating usage on a simple simulated dataset.

-   Parallelization of the "main loop" over the predictors in
    `covdepGE_main.R`). This is complicated by the `C++` code, however,
    two potential solutions are:

    -   [[StackOverflow
        suggestion](https://stackoverflow.com/questions/69789634/parallelization-of-rcpp-without-inline-creating-a-local-package?noredirect=1#comment123649680_69789634)]{style="color: blu"}

    -   [[RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)]{style="color: blu"}

    This is a finishing touch and **likely will not be implemented in
    the package by the end of the semester.**

::: thebibliography
1 Dasgupta, Sutanoy, et al. "An approximate Bayesian approach to
covariate dependent graphical modeling.\" 2021
:::
