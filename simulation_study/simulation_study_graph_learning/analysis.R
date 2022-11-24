rm(list = ls())
library(covdepGE)
library(loggle)
library(mgm)
library(varbvs)

load("disc_cov_dep_ntrials3_p11_n1_50_n2_50_lambda15_20221123_155056.Rda")
#load("disc_cov_dep_ntrials3_p11_n1_50_n2_50_lambda15_20221123_161228.Rda")

res <- results[setdiff(names(results), "sample_data")]
res <- lapply(res, lapply, `[`, c("sens", "spec"))
sens <- sapply(res, sapply, `[[`, "sens")
spec <- sapply(res, sapply, `[[`, "spec")
rowMeans(sens)
rowMeans(spec)
apply(sens, 1, sd)
apply(spec, 1, sd)
