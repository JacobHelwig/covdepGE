library(covdepGE)
library(JGL)
library(mgm)
# p = 25, n = 150
# load("/home/jacob.a.helwig/covdepGE/simulation_study/res_p25_n150_20220821_004616.Rda")

# p = 50, n = 150
# load("/home/jacob.a.helwig/covdepGE/simulation_study/res_p50_n150_20220821_004636.Rda")

# p = 100, n = 300
# load("/home/jacob.a.helwig/covdepGE/simulation_study/res_p100_n300_20220821_004700.Rda")

# p = 100, n = 600
# load("/home/jacob.a.helwig/covdepGE/simulation_study/res_p100_n600_20220821_004728.Rda")

# extract models
covdepGE_mods <- lapply(results, `[[`, "covdepGE")
jgl_mods <- lapply(results, `[[`, "JGL")
mgm_mods <- lapply(results, `[[`, "mgm")
mods <- list(covdepGE = covdepGE_mods,
             JGL = jgl_mods,
             mgm = mgm_mods)

# function for extracting
extractor <- function(lst, name) unlist(sapply(lst, `[[`, name))

# extract times, sensitivity, and specificity
times <- lapply(mods, extractor, "time")
sens <- lapply(mods, extractor, "sens")
spec <- lapply(mods, extractor, "spec")

# analyze
lapply(times, summary)
lapply(sens, summary)
lapply(spec, summary)
length(spec$covdepGE)
