rm(list = ls())
library(covdepGE)
library(loggle)
library(HeteroGGM)
library(mclust)

# initialize storage for results, time, and progress tracking
set.seed(1)
n_trials <- 100
results <- vector("list", n_trials + 1)
names(results) <- c(paste0("trial", 1:n_trials), "time")
results$time <- Sys.time()
pb <- txtProgressBar(0, n_trials, style = 3)

# function for surpressing cat
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


for (j in 1:n_trials){

  # generate the data and create storage for the models
  data <- generateData()
  n <- nrow(data$X)
  p <- ncol(data$X)
  trial <- vector("list", 4)
  names(trial) <- c("data", "covdepGE", "loggle", "HeteroGGM")
  prec_arr_dim <- c(p, rev(dim(data$X)))
  trial$data <- data

  # convert the true precision to an array and then to a graph; mask diagonal
  prec <- array(unlist(data$true_precision), prec_arr_dim)
  graph <- (prec != 0) * 1 + replicate(n, diag(rep(NA, p)) * 1)

  # get true number of edges and non-edges
  num_edge <- sum(graph, na.rm = T)
  num_non <- sum(graph == 0, na.rm = T)

  # fit each method, save details about results, time, etc.

  # covdepGE
  out_covdepGE <- suppressWarnings(suppressMessages(
    covdepGE(data$X, data$Z, parallel = T, num_workers = 5)))
  out_covdepGE$time <- as.numeric(out_covdepGE$model_details$elapsed, units = "secs")
  out_covdepGE$mem <- object.size(out_covdepGE)
  out_covdepGE$prec <- array(unlist(out_covdepGE$graphs$graphs), dim = prec_arr_dim)
  out_covdepGE$true_edge <- sum(out_covdepGE$prec == graph & graph == 1, na.rm = T)
  out_covdepGE$true_non <- sum(out_covdepGE$prec == graph & graph == 0, na.rm = T)
  out_covdepGE$sens <- out_covdepGE$true_edge / num_edge
  out_covdepGE$spec <- out_covdepGE$true_non / num_non
  trial$covdepGE <- out_covdepGE
  rm(list = "out_covdepGE")
  gc()

  # loggle
  num_workers <- parallel::detectCores()
  start <- Sys.time()
  out_loggle <- quiet(loggle.cv(t(data$X), num.thread = num_workers - 2))
  out_loggle$mem <- object.size(out_loggle)
  out_loggle <- out_loggle$cv.select.result
  out_loggle$time <- as.numeric(Sys.time() - start, units = "secs")
  out_loggle$prec <- array(unlist(lapply(lapply(
    out_loggle$adj.mat.opt, `-`, diag(p)), as.matrix)), dim = prec_arr_dim)
  out_loggle$true_edge <- sum(out_loggle$prec == graph & graph == 1, na.rm = T)
  out_loggle$true_non <- sum(out_loggle$prec == graph & graph == 0, na.rm = T)
  out_loggle$sens <-  out_loggle$true_edge / num_edge
  out_loggle$spec <-  out_loggle$true_non / num_non
  trial$loggle <- out_loggle
  rm(list = "out_loggle")
  gc()

  # HeteroGGM
  clust <- Mclust(data$Z, verbose = F)
  lambda <- genelambda.obo(lambda1_min = 0.01, lambda2_min = 0.2,
                           lambda3_min = 0.01)
  start <- Sys.time()
  out_hetGGM <- GGMPF(lambda, data$X + clust$classification * 10, clust$G)
  out_hetGGM$time <- as.numeric(Sys.time() - start, units = "secs")
  out_hetGGM$prec <- (out_hetGGM$Theta_hat.list[[out_hetGGM$Opt_num]] != 0) * 1
  out_hetGGM$prec <- out_hetGGM$prec - replicate(dim(out_hetGGM$prec)[3],
                                                 diag(p))
  out_hetGGM$prec <- out_hetGGM$prec[ , , out_hetGGM$member.list[[out_hetGGM$Opt_num]]]
  out_hetGGM$true_edge <- sum(out_hetGGM$prec == graph & graph == 1, na.rm = T)
  out_hetGGM$true_non <- sum(out_hetGGM$prec == graph & graph == 0, na.rm = T)
  out_hetGGM$sens <- out_hetGGM$true_edge / num_edge
  out_hetGGM$spec <- out_hetGGM$true_non / num_non
  trial$HeteroGGM <- out_hetGGM
  rm(list = "out_hetGGM")
  gc()

  # save the trial and update the progress bar
  results[[j]] <- trial
  setTxtProgressBar(pb, j)
}

# save the final time and the results
results$time <- Sys.time() - results$time
save(results, file = "competitor_res100.Rda")


## Analysis
load("dev/competitor_res100.Rda")

# extract each of the models
res <- results[setdiff(names(results), "time")]
covdepGE_models <- lapply(res, `[[`, "covdepGE")
loggle_models <- lapply(res, `[[`, "loggle")
hetGGM_models <- lapply(res, `[[`, "HeteroGGM")

# extract the sensitivties, specificities, and time
covdepGE_sens <- sapply(covdepGE_models, `[[`, "sens")
covdepGE_spec <- sapply(covdepGE_models, `[[`, "spec")
covdepGE_times <- sapply(covdepGE_models, `[[`, "time")
loggle_sens <- sapply(loggle_models, `[[`, "sens")
loggle_spec <- sapply(loggle_models, `[[`, "spec")
loggle_times <- sapply(loggle_models, `[[`, "time")
hetGGM_sens <- sapply(hetGGM_models, `[[`, "sens")
hetGGM_spec <- sapply(hetGGM_models, `[[`, "spec")
hetGGM_times <- sapply(hetGGM_models, `[[`, "time")

# sensitivity analysis
summary(covdepGE_sens)
summary(loggle_sens)
summary(hetGGM_sens)

summary(covdepGE_sens - loggle_sens)
summary(covdepGE_sens - hetGGM_sens)

# specificity analysis
summary(covdepGE_spec)
summary(loggle_spec)
summary(hetGGM_spec)

summary(covdepGE_spec - loggle_spec)
summary(covdepGE_spec - hetGGM_spec)

# time analysis
summary(covdepGE_times)
summary(loggle_times)
summary(hetGGM_times)

# put everything into a table
library(kableExtra)
perf <- apply(cbind(covdepGE_sens, loggle_sens, hetGGM_sens,
                    covdepGE_spec, loggle_spec, hetGGM_spec,
                    covdepGE_times, loggle_times, hetGGM_times),
              2, function(x) round(c(Mean = mean(x), SD = sd(x)), 3))
labs <- c("covdepGE", "loggle", "HetGGM")
perf_df <- data.frame(perf)
colnames(perf_df) <- rep(labs, 3)
kbl(perf_df, format = "latex", booktabs = T) %>%
  add_header_above(c("", "Sensitivity" = 3, "Specificity" = 3, "Time(s)" = 3)) %>%
  kable_styling(position = "center")

# plot
library(ggplot2)
library(ggsci)
# sensitivity
sens <- rbind.data.frame(
  cbind.data.frame(Sensitivity = covdepGE_sens, Package = "covdepGE"),
  cbind.data.frame(Sensitivity = loggle_sens, Package = "loggle"),
  cbind.data.frame(Sensitivity = hetGGM_sens, Package = "HeteroGGM"))
sens_meds <- aggregate(Sensitivity~Package, sens, function(x) round(median(x), 2))
ggplot(sens, aes(Package, Sensitivity, fill = Package)) + geom_boxplot(color = "#a7a7a7", position = position_dodge(1.2)) +
  coord_cartesian(ylim = c(0.7, 1)) +
  geom_text(data = sens_meds, aes(label = Sensitivity), color = "#d1d1d1", vjust = 1) +
  theme_classic() + scale_fill_manual(values = c("#500000", "#003C71", "#5B6236"))


p <- ggplot(data=mtcars, aes(x=factor(vs), y=mpg, fill=factor(am))) +
  geom_boxplot(position = position_dodge(width=0.9))

a <- aggregate(mpg ~ vs + am , mtcars, function(i) round(mean(i)))

p +  geom_text(data = a, aes(label = mpg))
