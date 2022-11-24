# USAGE: R CMD BATCH --no-save --no-restore "--args experiment='disc_cov_dep' p=11 n1=50 n2=50 lambda=15 n_trials=20" main.R main_all.Rout
# USAGE: R CMD BATCH --no-save --no-restore "--args experiment='cont_cov_dep' p=11 n1=50 n2=50 n3=50 n_trials=5" main.R cont_test.Rout
rm(list = ls())
source("models.R")
source("data.R")
source("utils.R")

# parse command line arguments
args <- (commandArgs(TRUE))
print(args)

# DEBUGGING
if (interactive()){
  args <- c("save_dir=./experiments", "experiment='cont_cov_dep'", "p=5", "n1=50", "n2=50", "n3=50", "n_trials=3")
  args <- c("save_dir=./experiments", "experiment='disc_cov_dep'", "p=11", "n1=50", "n2=50", "lambda=15", "n_trials=3")
}

if (length(args) > 0){
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
}else{
  stop("No arguments passed")
}

# check passed arguments
print(args)
if (!("experiment" %in% ls())){
  stop("Missing experiment")
}

# ensure experiment is recognized
disc_exps <- c("disc_cov_dep",
               "disc_cov_indep",
               "disc_cov_free")
cont_exps <- c("cont_cov_dep",
               "cont_multi_cov_dep")
experiment_choices <- c(disc_exps, cont_exps)
if (!(experiment %in% experiment_choices)){
  stop("Invalid experiment")
}

# check remaining args
expected_args <- c("save_dir", "p", "n1", "n2", "n_trials")
disc <- experiment %in% disc_exps
if (disc){
  expected_args <- c(expected_args, "lambda")
}else{
  expected_args <- c(expected_args, "n3")
}
if (!all(expected_args %in% ls())){
  stop("Missing args")
}

# create filename for saved results
set.seed(1)
(now <- format(Sys.time(), "%Y%m%d_%H%M%S"))
if (disc){
  file_suffix <- paste0("_lambda", lambda)
}else{
  file_suffix <- paste0("_n3_", n3)
}
filename <- paste0(experiment, "_ntrials", n_trials, "_p", p, "_n1_", n1,
                   "_n2_", n2, file_suffix, "_", now, ".Rda")
(filename <- file.path(save_dir, filename))

# create list for storing results
trial_list <- list(loggle = NA, mgm = NA, covdepGE = NA)
if (disc){
  trial_list$varbvs <- NA
}
results <- replicate(n_trials, trial_list, simplify = F)
names(results) <- c(paste0("trial", 1:n_trials))

# generate the data
set.seed(1)
if (experiment == "disc_cov_dep"){
  data_list <- replicate(n_trials, disc_cov_dep_data(p, n1, n2, lambda), F)
}else if (experiment == "cont_cov_dep"){
  data_list <- replicate(n_trials, cont_cov_dep_data(p, n1, n2, n3), F)
}else{
  stop("Experiment not yet implemented")
}

# perform trials
trials(data_list = data_list,
       results = results,
       filename = filename)
