
Here we demonstrate an application of `covdepGE` to real-world data
similar to that which was used in (1). In both cases, the protein
expression values for genes from patients with Breast Invasive Carcinoma
(BRCA) are modeled, with data sourced from The Cancer Genome Atlas
(TCGA). (1) considers the following genes :

- CTNNB1
- BRCA2
- MET
- **E-cadherin**
- **N-cadherin**
- **NFkB1**
- snail
- **STAT3**

Here, we consider the genes:

- STAT3_pY705
- E-Cadherin
- N-Cadherin
- NF-kB-p65_pS536

While (1) models the precision matrix for these genes as a function of
the FOXC2 gene, here we use FOXO3a. After fitting the model, we
visualize each of the estimated graphs (left column) alongside the
corresponding FOXO3a values (right column).

``` r
# BiocManager::install("RTCGA.BRCA")
library(covdepGE)
library(ggpubr)

# load data
dat <- RTCGA.RPPA::BRCA.RPPA
zvar <- "FOXO3a"
Z <- dat[ , zvar]
vars <- c("STAT3_pY705",
          "E-Cadherin",
          "N-Cadherin",
          "NF-kB-p65_pS536")
X <- dat[ , vars]

# fit model
out <- covdepGE(X, Z, parallel = T, num_workers = 4)
out
```

    ##                       Covariate Dependent Graphical Model
    ## 
    ## ELBO: -427829.27                                             # Unique Graphs: 5
    ## n: 410, variables: 4                       Hyperparameter grid size: 125 points
    ## Model fit completed in 5.767 secs

``` r
# returns graph visualization and histogram
plot_brca <- function(graph, z, xlim, ylim, color, labs, labz, bw){
  graph_viz <- matViz(graph, color2 = color) + 
    theme(legend.position='none', axis.text.x = element_text(angle=30,hjust=1,vjust=1.0)) + 
    scale_y_continuous(labels=labs, breaks=1:nrow(graph)) + 
    scale_x_continuous(labels=labs, breaks=1:nrow(graph))
  hist <- ggplot() + geom_histogram(aes(z), color='black', fill = color, binwidth = bw) + 
    xlab(labz) + theme_bw() + xlim(xlim) + ylim(c(0, ylim))
  list(graph_viz, hist)
}

# get graphs and corresponding observations
graphs <- lapply(out$graphs$unique_graphs, `[[`, "graph")
inds <- lapply(out$graphs$unique_graphs, `[[`, "indices")

# sort graphs by minimum Z value
ord <- order(sapply(inds, function(i) min(Z[i])))
colors <- ggsci::pal_nejm('default')(length(ord))

# plot results
get_graphs <- function(i){
  plot_brca(graph=graphs[[i]], z=Z[inds[[i]]], xlim=range(Z), ylim=25, color=colors[i], 
            labs=vars, labz=zvar, bw=0.04)
}
plots_list <- Reduce(c, lapply(ord, get_graphs))
ggarrange(plotlist = plots_list, ncol=2, nrow=length(ord))
```

![](TCGA_analysis_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Bibliography

1)  Dasgupta, Sutanoy, Peng Zhao, Jacob Helwig, Prasenjit Ghosh, Debdeep
    Pati, and Bani K. Mallick. “An Approximate Bayesian Approach to
    Covariate-dependent Graphical Modeling.” arXiv preprint
    arXiv:2303.08979 (2023).
