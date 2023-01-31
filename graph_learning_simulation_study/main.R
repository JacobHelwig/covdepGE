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
  # args <- c("save_dir='./experiments'", "experiment='cont_cov_dep'", "p=5", "n=50", "n2=50", "n3=50", "n_trials=3")
  # args <- c("save_dir='./experiments'", "experiment='cont_multi_cov_dep'", "p=11", "n=25", "n_trials=3")
  # args <- c("save_dir='./experiments'", "experiment='disc_cov_dep'", "p=11", "n1=80", "n2=20", "lambda=3", "n_trials=3")
  # args <- c("save_dir='./experiments'", "experiment='disc_cov_indep'", "p=11", "n1=50", "n2=50", "lambda=15", "n_trials=3")
  # args <- c("save_dir='./experiments'", "experiment='disc_cov_free'", "p=11", "n1=50", "n2=50", "lambda=15", "n_trials=3")
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
disc <- experiment %in% disc_exps
expected_args <- c("save_dir", "p", "n_trials")
if (experiment == "cont_multi_cov_dep"){
  expected_args <- c(expected_args, "n")
}else{
  expected_args <- c(expected_args, "n1", "n2")
  if (disc){
    expected_args <- c(expected_args, "lambda")
  }else{
    expected_args <- c(expected_args, "n3")
  }
}
if (!all(expected_args %in% ls())){
  stop("Missing args")
}

# check if any skips have been specified
if (!("skips" %in% ls())){
  skips <- NULL
}else{
  print(paste0(c("Skipping", skips), collapse = " "))
}

# check if a filename has been passed
if ("filename" %in% ls()){

  # a file has been specified; load the associated results
  load(filename)

}else{

  # a filename has not been specified

  # create list for storing results
  trial_list <- list(mgm = NA, covdepGE = NA)
  if (disc){
    trial_list$varbvs <- NA
    if (experiment == "disc_cov_free"){
      trial_list <- list(varbvs = NA, covdepGE = NA)
    }
  }else{
    trial_list$loggle <- NA
    if (experiment == "cont_multi_cov_dep"){
      trial_list$covdepGE_sortZ <- NA
    }
  }
  results <- replicate(n_trials, trial_list, simplify = F)
  names(results) <- c(paste0("trial", 1:n_trials))

  # create filename for saved results
  (now <- format(Sys.time(), "%Y%m%d_%H%M%S"))
  filename <- paste0(c(setdiff(names(trial_list), skips), ""), collapse = "_")
  filename <- paste0(filename, experiment, "_ntrials", n_trials, "_p", p)
  if (experiment == "cont_multi_cov_dep"){
    filename <- paste0(filename, "_n_", n, "_", now, ".Rda")
  }else{
    if (disc){
      file_suffix <- paste0("_lambda", lambda)
    }else{
      file_suffix <- paste0("_n3_", n3)
    }
    filename <- paste0(filename, "_n1_", n1, "_n2_", n2, file_suffix, "_", now,
                       ".Rda")
  }
  (filename <- file.path(save_dir, filename))
}

# generate the data
set.seed(1)
if (experiment == "disc_cov_dep"){
  data_list <- replicate(n_trials, disc_cov_dep_data(p, n1, n2, lambda), F)
}else if (experiment %in% c("disc_cov_indep", "disc_cov_free")){
  data_list <- replicate(n_trials,
                         disc_cov_dep_data(p, n1, n2, lambda, independent = T), F)
  if (experiment == "disc_cov_free"){
    for (j in 1:n_trials){
      data_list[[j]]$Z <- rep(1, n1 + n2)
    }
  }
}else if (experiment == "cont_cov_dep"){
  data_list <- replicate(n_trials, cont_cov_dep_data(p, n1, n2, n3), F)
}else{
  data_list <- replicate(n_trials, cont_multi_cov_dep_data(p, n), F)
}

# DEBUGGING
if (interactive()){
  X_ <- data_list[[1]]$X
  Z_ <- data_list[[1]]$Z
  true_ <- data_list[[1]]$true_precision
  out_ <- covdepGE(X_, Z_)
  out_$str <- array(unlist(out_$graphs$graphs), dim = c(p, p, nrow(X_)))
  (perf_ <- eval_est(out_$str, true_))
}

# perform trials
trials(data_list = data_list,
       results = results,
       filename = filename,
       skips = skips)
