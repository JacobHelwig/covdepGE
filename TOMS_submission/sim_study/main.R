# R CMD BATCH --no-save --no-restore "--args save_dir='./experiments/z2' experiment='cont_multi_cov_dep' p=100 n=25 n_trials=50 skips=c('JGL','mgm','covdepGE') trial_skips=1:49" main.R ./experiments/z2/cont_multi_cov_dep_ntrials50_p100_n_25_2023-02-04_13-15.Rout
rm(list = ls())
source("models.R")
source("data.R")
source("utils.R")

# parse command line arguments
args <- (commandArgs(TRUE))
print(args)

# DEBUGGING
if (interactive()){
  # args <- c("save_dir='./experiments/z2'", "experiment='cont_multi_cov_dep'", "p=10", "n=10", "n_trials=10")
  # args <- c("save_dir='./experiments'", "experiment='cont_cov_dep'", "p=3", "n1=5", "n2=5", "n3=5", "n_trials=1")
  # args <- c("save_dir='./experiments'", "experiment='cont_multi_cov_dep'", "p=3", "n=2", "n_trials=1")
  # args <- c("save_dir='./experiments'", "experiment='cont_cov_dep_sine'", "p=6", "n1=75", "n2=75", "n3=75", "n_trials=1")
  args <- c("save_dir='./experiments'", "experiment='cont_4_cov_dep_data'", "p=10", "n=225", "n_trials=1")
  # args <- c("save_dir='./experiments'", "experiment='cont_cov_dep'", "p=5", "n1=50", "n2=50", "n3=50", "n_trials=1")
  # args <- c("save_dir='./experiments'", "experiment='cont_multi_cov_dep'", "p=11", "n=25", "n_trials=1")
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
experiment_choices <- c("cont_cov_dep", "cont_cov_dep_sine", "cont_multi_cov_dep", "cont_4_cov_dep_data")
if (!(experiment %in% experiment_choices)){
  stop("Invalid experiment")
}

# check remaining args
expected_args <- c("save_dir", "p", "n_trials")
if (experiment %in% c("cont_multi_cov_dep", "cont_4_cov_dep_data")){
  expected_args <- c(expected_args, "n")
}else{
  expected_args <- c(expected_args, "n1", "n2", "n3")
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
if (!("trial_skips" %in% ls())){
  trial_skips <- NULL
}else{
  print(paste0(c("Skipping trials", trial_skips), collapse = " "))
}
# check if a HP method or max_iter_grid for covdepGE has been specified
if (!("hp_method" %in% ls())) hp_method <- "hybrid"
if (!("max_iter_grid" %in% ls())) max_iter_grid <- 10
print(paste0(c("covdepGE HP specification scheme: ", hp_method), collapse = " "))

# check if a filename has been passed
if ("filename" %in% ls()){

  # a file has been specified; load the associated results
  load(filename)

}else{

  # a filename has not been specified

  # create list for storing results
  trial_list <- list(covdepGE = NA, JGL = NA, mgm = NA)
  if (experiment == "cont_multi_cov_dep"){
    trial_list$covdepGE_sortZ <- NA
  }

  results <- replicate(n_trials, trial_list, simplify = F)
  names(results) <- c(paste0("trial", 1:n_trials))

  # create filename for saved results
  (now <- format(Sys.time(), "%Y%m%d_%H%M%S"))
  filename <- paste0(c(setdiff(names(trial_list), skips), ""), collapse = "_")
  filename <- gsub("covdepGE" , paste0("covdepGE", hp_method) ,filename)
  filename <- paste0(filename, experiment, "_ntrials", n_trials, "_p", p)
  if (experiment %in% c("cont_multi_cov_dep", "cont_4_cov_dep_data")){
    filename <- paste0(filename, "_n_", n, "_", now, ".Rda")
  }else{
    filename <- paste0(filename, "_n1_", n1, "_n2_", n2, "_n3_", n3, "_", now,
                       ".Rda")
  }
  (filename <- file.path(save_dir, filename))
}

# generate the data
set.seed(1)
set.seed(2)
if (experiment == "cont_cov_dep"){
  data_list <- replicate(n_trials, cont_cov_dep_data(p, n1, n2, n3), F)
}else if (experiment == "cont_cov_dep_sine"){
  data_list <- replicate(n_trials, cont_cov_dep_sine_data(p, n1, n2, n3), F)
}else if (experiment == "cont_4_cov_dep_data"){
  data_list <- replicate(n_trials, cont_4_cov_dep_data(p, n), F)
}else{
  data_list <- replicate(n_trials, cont_multi_cov_dep_data(p, n), F)
}

# DEBUGGING
if (interactive()){
  X_ <- data_list[[1]]$X
  Z_ <- data_list[[1]]$Z
  true_ <- data_list[[1]]$true_precision
  graph_viz(lapply(true_, `!=`, 0))
  out_ <- covdepGE(X_, Z_)
  plot(out_)
  out_$str <- array(unlist(out_$graphs$graphs), dim = c(p, p, nrow(X_)))
  (perf_ <- eval_est(out_$str, true_))
}

# perform trials
trials(data_list = data_list,
       results = results,
       filename = filename,
       skips = skips,
       trial_skips = trial_skips,
       hp_method = hp_method,
       max_iter_grid = max_iter_grid)

