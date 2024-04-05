# ------------------------------------------------------------------------------
# Baseline comparisons
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("data.R")
source("utils.R")
library(kableExtra)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(latex2exp)
library(scales)
windowsFonts("Times" = windowsFont("Times"))
windowsFonts("Courier" = windowsFont("Courier New"))

# init list for recording statistics for each graph; only for q=1
graphstats_list <- vector("list", 4)
names(graphstats_list) <- c("10", "25", "50", "100")
graphstats_list <- list(nl=graphstats_list, pwl=graphstats_list)

# init list for counting number of graphs estimated in each experiment
graph_cts <- list(
  cont_cov_dep_=list(`10`=rep(3,50), `25`=rep(3,50), `50`=rep(3,50), `100`=rep(3,50)),
  cont_cov_dep_sine_=list(`10`=rep(3,50), `25`=rep(3,50), `50`=rep(3,50), `100`=rep(3,50)),
  cont_multi_cov_dep_=list(`10`=rep(4,50), `25`=rep(4,50), `50`=rep(4,50), `100`=rep(4,50)),
  cont_4_cov_dep_=list(`10`=c(16, 15, 14, 16, 16, 16, 16, 15, 15, 15, 14, 14, 16, 16, 16,
                              15, 16, 15, 16, 16, 16, 15, 15, 16, 14, 16, 15, 16, 16, 16, 15,
                              16, 15, 15, 16, 16, 16, 13, 16, 16, 14, 16, 16, 16, 14, 16, 15,
                              14, 15, 15),
                       `25`=c(16, 15, 16, 15, 15, 16, 16, 16, 16, 16, 16, 15, 16, 16, 15,
                              16, 16, 15, 15, 16, 16, 16, 16, 16, 16, 15, 16, 15, 15, 15, 15,
                              16, 14, 16, 16, 16, 16, 16, 15, 16, 16, 16, 16, 16, 16, 16, 16,
                              16, 16, 16),
                       `50`=c(16, 16, 16, 16, 16, 16, 15, 15, 16, 16, 16, 15, 16, 16, 15,
                              16, 14, 15, 14, 16, 15, 16, 16, 16, 16, 14, 16, 16, 16, 14, 16,
                              14, 16, 15, 16, 16, 15, 16, 15, 16, 15, 14, 15, 14, 16, 16, 15,
                              15, 14, 16),
                       `100`=c(16, 16, 15, 16, 16, 15, 15, 16, 16, 16, 16, 15, 15, 15, 16,
                               16, 16, 15, 16, 16, 14, 16, 15, 15, 16, 16, 16, 16, 16, 16, 15,
                               16, 15, 16, 15, 16, 15, 15, 16, 16, 16, 16, 16, 16, 15, 15, 16,
                               16, 16, 16))
)
ngraphs <- lapply(graph_cts, lapply, function(x) NULL)

# init list for recording number of clusters estimated by mclust
subgroups_list <- list()

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each

# results config; med outputs median in place of mean and univariate is for
# switching setting from q=1->q=2
med <- F # median aggregate results
univariate <- T # T for q=1, F for q=2,4
sine <- T # q=1 w non-linear covariate
four <- F # q=4
seq <- F # aggregate times for sequentials
subgr <- T # analyze performance for each graph for comparing sensitivity/specificity in each setting
if (univariate){

  # univariate extraneous covariate
  exper <- "cont_cov_dep_"
  methods <- c("covdepGE" , "JGL", "loggle", "mgm")
  n <- 75
  samp_str <- paste0("_n1_", n, "_n2_", n, "_n3_", n)
  path <- "./experiments/z1"
  if (sine){
    path <- paste0(path, "_sine")
    exper <- paste0(exper, "sine_")
  }else if (seq){
    seq_path <- paste0(path, "_seq")
    seq_files <- list.files(seq_path)
    methods <- c("covdepGE", "covdepGE_seq")
  }
}else{

  # multivariate extraneous covariate
  exper <- ifelse(four, "cont_4_cov_dep_", "cont_multi_cov_dep_")
  methods <- c("covdepGE" , "covdepGE_sortZ", "JGL", "loggle", "mgm")
  n <- 25 + 200 * four
  samp_str <- paste0("_n_", n)
  path <- ifelse(four, "./experiments/z4", "./experiments/z2")
}

dims <- c(10, 25, 50, 100)
ntrials <- 50
trial_str <- paste0("ntrials", ntrials, "_")
files <- list.files(path)
loggle_path <- paste0(path, '/loggle')
loggle_files <- list.files(loggle_path)
prec <- 2
df <- subgroups <- htest <- vector("list", length(dims))
names(df) <- names(subgroups) <- names(htest) <- dims
df_num <- list(sens=data.frame(), spec=data.frame(), time=data.frame())

if (subgr){
  stats_name <- ifelse(sine, "nl", "pwl")

  set.seed(1)
  sine_data <- cont_cov_dep_sine_data(p=3, n1=75, n2=75, n3=75)
  table(cut(c(sine_data$Z), c(-3, -2, -1, 1, 2, 3), right=F)) # TODO: update description in paper

  # init list for recording statistics for specific subgraphs; only for q=1 settings
  subg <- list(pwl=list(), nl=list())
  get_lab <- function(g) paste0("CDS", which(sapply(ugraphs, function(gr) isTRUE(all.equal(g, gr)))))
  for (dim in as.character(dims)){

    # get precision
    subg$pwl[[dim]] <- list(true=cont_cov_dep_data(p=as.numeric(dim), n1=75, n2=75, n3=75)$true_precision)
    subg$nl[[dim]] <- list(true=cont_cov_dep_sine_data(p=as.numeric(dim), n1=75, n2=75, n3=75)$true_precision)

    # get graphs
    graphs_pwl <- lapply(lapply(subg$pwl[[dim]]$true, `!=`, 0), `-`, diag(dim))
    ugraphs <- unique(graphs_pwl)
    graphs_nl <- lapply(lapply(subg$nl[[dim]]$true, `!=`, 0), `-`, diag(dim))

    # get labels (CDS1,2,3)
    subg$pwl[[dim]]$label <- sapply(graphs_pwl, get_lab)
    subg$nl[[dim]]$label <- sapply(graphs_nl, get_lab)
  }
}
loggle <- F
p <- "10"

for (p in as.character(dims)){

  # format file name
  exp_name <- paste0(exper, trial_str, "p", p, samp_str)

  # load loggle/seq results
  if (!seq){
    file_name <- loggle_files[grepl(exp_name, loggle_files) & endsWith(loggle_files, ".Rda")]
    if (length(file_name) != 1) stop(paste(length(file_name), "files found for", exp_name))
    file_path <- file.path(loggle_path, file_name)
    load(file_path)
    loggle_results <- results
  }else{
    file_name <- seq_files[grepl(exp_name, seq_files) & endsWith(seq_files, ".Rda")]
    if (length(file_name) != 1) stop(paste(length(file_name), "files found for", exp_name))
    file_path <- file.path(seq_path, file_name)
    load(file_path)
    seq_results <- results
  }

  # load results
  # mean(sapply(results, function(trial)length(trial$covdepGE$graphs$unique_graphs)))
  if (!loggle){
    file_name <- files[grepl(exp_name, files) & endsWith(files, ".Rda")]
    if (length(file_name) != 1) stop(paste(length(file_name), "files found for", exp_name))
    file_path <- file.path(path, file_name)
    load(file_path)
  }
  print(paste0("n: ", results$sample_data[1], ", p: ", results$sample_data[2]))
  results <- results[setdiff(names(results), "sample_data")]

  # add loggle/seq results
  if (!seq){
    for (i in 1:length(results)){
      results[[i]]$loggle <- loggle_results[[i]]$loggle
    }
  }else{
    for (i in 1:length(results)){
      results[[i]]$covdepGE_seq <- seq_results[[i]]$covdepGE_seq
    }
  }

  # get graph stats
  if (subgr){

    graphs <- vector("list", length(results))

    # CDS masks
    true_graphs <- subg[[stats_name]][[p]]$true
    labels <- subg[[stats_name]][[p]]$label
    ints <- c("CDS1", "CDS2", "CDS3")
    masks <- lapply(ints, `==`, labels)
    names(masks) <- ints

    # reconstruct graphs and get stats
    int_res <- vector("list", 3)
    names(int_res) <- sapply(1:3, function(i) paste0("CDS", i))
    stats_list <- replicate(50, int_res, F)
    for (i in 1:length(results)){
      graphs[[i]] <- vector("list", 225)
      preds <- results[[i]]$covdepGE$graphs$unique_graphs
      for (pred in preds){
        graphs[[i]][pred$indices] <- replicate(length(pred$indices), as.matrix(pred$graph), F)
      }

      graphs[[i]] <- array(unlist(graphs[[i]]), c(p, p, 3 * n))

      # check for correct reconstruction
      full_stats <- eval_est(graphs[[i]], true_graphs)
      if (!all.equal(unlist(results[[i]]$covdepGE[names(full_stats)]), unlist(full_stats))){
        stop(paste0(i, " not equal"))
      }

      # get stats in each interval
      for (int in ints){
        mask <- masks[[int]]
        stats_list[[i]][[int]] <- eval_est(graphs[[i]][,,mask], true_graphs[mask])
      }

      # check for construct calculation
      if (stats_name == "pwl"){
        full_stats2 <- as.list(rowMeans(matrix(unlist(stats_list[[i]]), ncol=3)))
        names(full_stats2) <- names(full_stats)
        full_stats <- unlist(full_stats[setdiff(names(full_stats), c("sens", "spec"))])
        full_stats2 <- unlist(full_stats2[setdiff(names(full_stats2), c("sens", "spec"))])
        if (!all.equal(full_stats, full_stats2)){
          stop(paste0(i, " not equal"))
        }
      }

    }

    # process sensitivity results
    sens <- sapply(stats_list, sapply, `[[`, "sens")
    sens_mean <- rowMeans(sens) * 100
    sens_sd  <- apply(sens * 100, 1, sd)
    graphstats_list[[stats_name]][[p]]$sens <- list(data=sens, mean=sens_mean, sd=sens_sd)

    # process specificity results
    spec <- sapply(stats_list, sapply, `[[`, "spec")
    spec_mean <- rowMeans(spec) * 100
    spec_sd  <- apply(spec * 100, 1, sd)
    graphstats_list[[stats_name]][[p]]$spec <- list(data=spec, mean=spec_mean, sd=spec_sd)

  }

  # extract number of graphs
  ngraphs[[exper]][[p]] <- list(pred=sapply(results, function(trial)length(trial$covdepGE$graphs$unique_graphs)),
                                true=graph_cts[[exper]][[p]])
  ngraphs[[exper]][[p]]$perc <- ngraphs[[exper]][[p]]$pred / ngraphs[[exper]][[p]]$true

  # out <- results$trial2$covdepGE
  # for (j in 1:length(out$graphs$unique_graphs)){
  #   out$graphs$unique_graphs[[j]]$graph <- as.matrix(out$graphs$unique_graphs[[j]]$graph)
  # }
  # plot(out)

  # process sensitivity results
  sens <- sapply(results, sapply, `[[`, "sens") * 100
  sens_mean <- sens_mean_hyp <- rowMeans(sens)
  if (med) sens_mean <- apply(sens, 1, median)
  max_sens_ind <- which.max(sens_mean)
  sens_mean_str <- sprintf(paste0("%.", prec, "f"), sens_mean)
  sens_mean_str[max_sens_ind] <- paste0("\\mathbf{", sens_mean_str[max_sens_ind], "}")
  sens_sd  <- apply(sens, 1, sd)
  sens_se <- sens_sd / sqrt(ntrials)
  sens_se_str <- sprintf(paste0("%.", prec, "f"), sens_se)
  if (med) sens_se_str <- sprintf(paste0("%.", prec, "f"), apply(sens, 1, IQR))
  sens_str <- paste0(sens_mean_str, " (\\new{", sens_se_str, "})")

  # hypothesis testing
  if (!seq){
    max_baseline <- names(which.max(sens_mean_hyp[!sapply(names(sens_mean_hyp), function(x) grepl("covdepGE", x))]))
    t_test <- t.test(sens['covdepGE',], sens[max_baseline,], alternative="greater")
    pval <- sprintf(paste0("%.", prec + 1, "f"), t_test$p.value)
    pval <- ifelse(round(t_test$p.value,prec+1) == 0, paste0(c("<0.",rep(0, prec), 1), collapse=""), pval)
    # pval <- ifelse(t_test$p.value < 0.05, paste0("\\mathbf{", pval, "}"), pval)
    pval <- paste0("$", pval, "$")
    stat <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$statistic), "$")
    degf <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$parameter), "$")
    mean_covdepGE <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$estimate[1] * 100), "$")
    # mean_JGL <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$estimate[2] * 100), "$")
    se <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$stderr), "$")
    htest[[p]] <- data.frame(covdepGE=paste0("$", sens_str[row.names(sens) == "covdepGE"], "$"),
                             baseline_mean=paste0("$", sens_str[row.names(sens) == max_baseline], "$"),
                             baseline=paste0("\\texttt{", max_baseline, "}"),
                             stderr=se, stat=stat, df=degf, p=pval)
  }

  # process specificity results
  spec <- sapply(results, sapply, `[[`, "spec") * 100
  spec_mean <- rowMeans(spec)
  if (med) spec_mean <- apply(spec, 1, median)
  max_spec_ind <- which.max(spec_mean)
  spec_mean_str <- sprintf(paste0("%.", prec, "f"), spec_mean)
  spec_mean_str[max_spec_ind] <- paste0("\\mathbf{", spec_mean_str[max_spec_ind], "}")
  spec_sd  <- apply(spec, 1, sd)
  spec_se <- spec_sd / sqrt(ntrials)
  spec_se_str  <- sprintf(paste0("%.", prec, "f"), spec_se)
  if (med) spec_se_str <- sprintf(paste0("%.", prec, "f"), apply(spec, 1, IQR))
  spec_str <- paste0(spec_mean_str, " (\\new{", spec_se_str, "})")

  # process time results
  time <- sapply(results, sapply, `[[`, "time")
  time_mean <- rowMeans(time)
  if (seq){
    time_ratios <- time['covdepGE_seq',] / time['covdepGE',]
    ratio_mean <- sprintf(paste0("%.", prec, "f"), mean(time_ratios))
    ratio_sd <- sd(time_ratios)
    ratio_se <- ratio_sd / sqrt(ntrials)
    ratio_se <- sprintf(paste0("%.", prec, "f"), ratio_se)
    ratio_str <- paste0(ratio_mean, " (\\new{", ratio_se, "})")
  }
  if (med) time_mean <- apply(time, 1, median)
  min_time_ind <- which.min(time_mean)
  time_mean_str <- sprintf(paste0("%.", prec, "f"), time_mean)
  time_mean_str[min_time_ind] <- paste0("\\mathbf{", time_mean_str[min_time_ind], "}")
  time_sd <- apply(time, 1, sd)
  time_se <- time_sd / sqrt(ntrials)
  time_se_str  <- sprintf(paste0("%.", prec, "f"), time_se)
  if (med) time_se_str <- sprintf(paste0("%.", prec, "f"), apply(time, 1, IQR))
  time_str <- paste0(time_mean_str, " (\\new{", time_se_str, "})")

  # combine summary strings
  perf_str <- cbind(sens_str, spec_str, time_str)
  perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
  row.names(perf_str) <- row.names(spec)

  # combine summary numerics
  df_num[["sens"]] <- rbind.data.frame(df_num[["sens"]],data.frame(p=p, Method=names(sens_mean), mean=sens_mean, se=sens_se))
  df_num[["spec"]] <- rbind.data.frame(df_num[["spec"]],data.frame(p=p, Method=names(spec_mean), mean=spec_mean, se=spec_se))
  df_num[["time"]] <- rbind.data.frame(df_num[["time"]],data.frame(p=p, Method=names(time_mean), mean=time_mean, se=time_se))

  # create storage
  df_exp <- data.frame(p = p, method = methods, sens = NA, spec = NA, time = NA)
  df_exp[ , c("sens", "spec", "time")] <- perf_str[df_exp$method, ]
  if (seq){
    df_exp$time_ratio <- ratio_str
  }
  df[[p]] <- df_exp

  subgroups[[p]] <- sapply(results, function(trial) length(unique(trial$JGL$classification)))

  rm("results")
}
if (F){
  colors <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF") # https://nanx.me/ggsci/reference/pal_nejm.html
  exp_map <- list("$\\textit{q}=1$, PWL", "$\\textit{q}=1$, NL", "$\\textit{q}=2$", "$\\textit{q}=4$")
  names(exp_map) <- names(ngraphs)
  plots <- lapply(1:length(ngraphs), function(exp_ind) lapply(1:4, function(p_ind) ggplot() +
      geom_histogram(aes(x = ngraphs[[exp_ind]][[p_ind]]$pred),
                     color = "black", fill = colors[exp_ind], binwidth = ifelse(max(ngraphs[[exp_ind]][[p_ind]]$pred) < 10, 1, ifelse(max(ngraphs[[exp_ind]][[p_ind]]$pred) < 40, 2, 4))) +
      theme_pubclean() +
      theme(text = element_text(family = "Times", size = 14),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(TeX(paste0("$\\textit{p}=", names(ngraphs[[exp_ind]])[p_ind], "$, ", exp_map[[exp_ind]]))) +
      labs(x = TeX("Number of Unique Graphs")) + scale_y_continuous(breaks = scales::pretty_breaks()) + scale_x_continuous(breaks = scales::pretty_breaks())
  ))
  all_plots <- lapply(unlist(plots, recursive = F), function(g) g + rremove("ylab") + rremove("xlab"))
  # all_plots <- unlist(plots, recursive = F)
  arplots <- ggarrange(plotlist=all_plots, nrow=4, ncol=4)
  annotate_figure(arplots, left=text_grob("Number of Trials",  size = 18, family="Times", rot = 90),
                  bottom=text_grob("Number of Unique Graphs",  size = 18, family="Times"))
  ggsave("plots/unique_graphs.pdf", height = 10, width = 10)
}
if (subgr){
  graphstats_list_backup <- graphstats_list

  # https://stackoverflow.com/questions/18165863/multirow-axis-labels-with-nested-grouping-variables

  # # gather data for each graph into a dataframe where the columns are the
  # # covariate type, p, CDS index, and the data for that choice
  # stats <- lapply(graphstats_list, lapply, lapply, `[[`, "data")
  # sens <- lapply(stats, lapply, `[[`, "sens")
  # sens <- Reduce(rbind, lapply(c("nl", "pwl") , function(cov_type) Reduce(rbind, lapply(names(sens[[cov_type]]), function(p) data.frame(cov_type, p,t(sens[[cov_type]][[p]]))))))
  # sens <- reshape(sens, direction="long", varying=3:5, v.names="CDS")
  # names(sens) <-c("cov_type", "p", "CDS", "sens", "id")
  # sens <- sens[, setdiff(names(sens), "id")]
  # sens_pwl <- sens[sens$cov_type == "pwl",]
  # sens_nl <- sens[sens$cov_type == "nl",]
  # ggplot(data = sens, aes(x = p, y = sens, fill = cov_type)) +
  #   geom_boxplot() +
  #   facet_wrap(~CDS, strip.position = "bottom", scales = "free_x") +
  #   theme(panel.spacing = unit(0, "lines"),
  #         strip.background = element_blank(),
  #         strip.placement = "outside")
  #
  #
  # # gather data for each graph into a dataframe where the columns are the
  # # covariate type, p, CDS index, and the data for that choice
  # stats <- lapply(c("mean", "sd"), function(stat) lapply(graphstats_list, lapply, lapply, `[[`, stat))
  # names(stats) <- c("mean", "sd")
  # sens <- lapply(stats, lapply, lapply, `[[`, "sens")
  # sens <- lapply(sens, function(stat_list) Reduce(rbind, lapply(c("nl", "pwl") , function(cov_type) Reduce(rbind, lapply(names(stat_list[[cov_type]]), function(p) data.frame(cov_type, p,t(stat_list[[cov_type]][[p]])))))))
  # sens <- lapply(sens, function(stat_list) reshape(stat_list, direction="long", varying=3:5, v.names="CDS"))
  # names(sens$mean) <- names(sens$sd) <- c("cov_type", "p", "CDS", "sens", "id")
  # sens$mean <- sens$mean[, setdiff(names(sens$mean), "id")]
  # sens$sd <- sens$sd[, setdiff(names(sens$sd), "id")]
  # sens$mean$p <- as.numeric(sens$mean$p)
  #
  # ggplot(data = sens$mean, aes(x = CDS, y = sens, fill = cov_type)) +
  #   geom_bar(stat = "identity", width = 1, position = 'dodge') +
  #   facet_wrap(~p, strip.position = "bottom", scales = "free_x") +
  #   theme(panel.spacing = unit(0, "lines"),
  #         strip.background = element_blank(),
  #         strip.placement = "outside")

  # gather data for each graph into a dataframe where the columns are the
  # covariate type, p, CDS index, and the data for that choice
  stats <- lapply(c("mean", "sd"), function(stat) lapply(graphstats_list, lapply, lapply, `[[`, stat))
  names(stats) <- c("mean", "sd")
  sens <- lapply(stats, lapply, lapply, `[[`, "sens")
  sens <- lapply(sens, function(stat_list) Reduce(rbind, lapply(c("nl", "pwl") , function(cov_type) Reduce(rbind, lapply(names(stat_list[[cov_type]]), function(p) data.frame(cov_type, p,t(stat_list[[cov_type]][[p]])))))))
  sens <- lapply(sens, function(stat_list) reshape(stat_list, direction="long", varying=3:5, v.names="CDS"))
  names(sens$mean) <- c("cov_type", "p", "CDS", "mean", "id")
  names(sens$sd) <- c("cov_type", "p", "CDS", "sd", "id")
  sens$mean <- sens$mean[, setdiff(names(sens$mean), "id")]
  sens$sd <- sens$sd[, setdiff(names(sens$sd), "id")]
  n <- ncol(graphstats_list$nl$`10`$sens$data)
  sens <- merge(sens$mean, sens$sd)
  sens$se <- sens$sd / sqrt(n)
  sens$p <- factor(sens$p, levels = c(10, 25, 50, 100), labels=c("10"=parse(text=TeX("$\\textit{p}=10$")),
                                    "25"=parse(text=TeX("$\\textit{p}=25$")),
                                    "50"=parse(text=TeX("$\\textit{p}=50$")),
                                    "100"=parse(text=TeX("$\\textit{p}=100$"))))
  sens$CDS <- paste("CDS", sens$CDS)
  sens$cov_type <- factor(toupper(sens$cov_type), levels = c("PWL", "NL"))

  windowsFonts("Times" = windowsFont("Times"))
  library(ggsci)
  library(scales)
  ggplot(data = sens, aes(x = CDS, y = mean, fill = cov_type)) +
    geom_bar(stat = "identity", width = 1, position = 'dodge', color="black") +
    geom_errorbar(aes(ymin=mean-2 * se, ymax=mean+2 * se),
                  # size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(1)) +
    facet_wrap(~p, strip.position = "bottom", scales = "free_x", labeller=label_parsed) +
    theme_pubclean() +
    theme(text = element_text(family = "Times", size = 18),
          plot.title = element_text(hjust = 0.5), panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          strip.placement = "outside", legend.title=element_blank()) +
    ggtitle(TeX(paste0("Sensitivity by CDS"))) +
    scale_y_continuous(limits = c(40, 100), oob=rescale_none) +
    labs(x=NULL, y="Sensitivity (%)", legend=NULL) + ggsci::scale_fill_nejm()
  ggsave("plots/q1_sens_black_outl.pdf", height = 8, width = 12)

}

if (seq){
  df <- Reduce(rbind, df)
  df <- df[,c("p", "method", "time", "time_ratio")]
  df <- data.frame(p = df[df$method == "covdepGE", "p"],
                   ptime = df[df$method == "covdepGE","time"],
                   stime = df[df$method == "covdepGE_seq","time"],
                   ratio = df[df$method == "covdepGE_seq","time_ratio"])
  colnames(df) <- colnames(df) <- c("$p$", "Parallel Time (seconds)", "Sequential Time (seconds)", "Ratio")
  kbl(df, format = "latex", booktabs = T, escape = FALSE) # %>%
    collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")
}

htest_df <- cbind(p=names(htest), Reduce(rbind, htest))
colnames(htest_df) <- c("$p$", "\\texttt{covdepGE}", "Baseline", "Baseline Name", "Standard Error", "T-statistic", "DF", "p-value")
kbl(htest_df, format = "latex", booktabs = T, escape = FALSE)

for (p in names(df)){
  df[[p]][df[[p]]$method == "loggle", setdiff(names(df[[p]]), 'p')] <- sapply(
    df[[p]][df[[p]]$method == "loggle", setdiff(names(df[[p]]), 'p')], function(x)paste0("\\new{", x, "}"))
}
df <- Reduce(rbind, df)
df$method <- gsub("_sortZ", "\\\\_time", df$method)
df$method <- paste0("\\texttt{", df$method, "}")
colnames(df) <- c("$p$", "Method", "Sensitivity$(\\%)$", "Specificity$(\\%)$", "Time (seconds)")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

for (metric in names(df_num)){
  df_num[[metric]]$p <- factor(df_num[[metric]]$p, levels=c(10, 25, 50, 100))
  df_num[[metric]]$Method <- gsub("_sortZ", "_time", df_num[[metric]]$Method)
  # df_num[[metric]]$Method <- paste0("\\textit{", df_num[[metric]]$Method, "}")
  levels <- unique(df_num[[metric]]$Method)[order(unique(df_num[[metric]]$Method))]
  colors <- ggsci::pal_nejm()
  if (!univariate){
    pal_colors <- colors(5)[c(1,5,2,3,4)]
  }else{
    pal_colors <- colors(4)
  }
  df_num[[metric]]$Method <- factor(df_num[[metric]]$Method, levels = levels)
  dodge <- unique(as.numeric(df_num[[metric]]$Method)/20)-mean(unique(as.numeric(df_num[[metric]]$Method)/20))
  df_num[[metric]]$p2 <- as.numeric(df_num[[metric]]$p) + dodge
}

plt_res <- function(data, ylab, xlab="", log=F, breaks=NULL, sz=0.75){
  plt <- ggplot(data, aes(x=p2, y=mean, group=Method, color=Method)) +
    theme_pubclean() + scale_color_manual(values = pal_colors, drop=FALSE) +
    theme(plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          text = element_text(family = "Times", size=18),
          legend.title=element_blank(), legend.key=element_blank(),
          legend.text=element_text(family="Courier", size=18)) +
    labs(x=xlab, y=ylab) +
    lapply(rev(levels), function(method) {
      list(
        geom_errorbar(data = ~ subset(., Method == method),
                      aes(ymin = mean - ifelse((log & mean - 2 * se < 0), se, 2 * se), ymax = mean + 2 * se),
                      width = 0.2, size = sz),
        geom_point(data = ~ subset(., Method == method), size = 2, aes()),
        geom_line(data = ~ subset(., Method == method), size = sz)
      )
    }) + scale_x_continuous(labels=c(10, 25, 50, 100), breaks = 1:4)

  if (log){
    plt <- plt + scale_y_continuous(trans='log10')
  }
  if (!is.null(breaks)){
    plt <- plt + scale_y_continuous(breaks = breaks)
  }
  plt
}

sens_breaks <- NULL # ifelse(univariate, ifelse(sine, NULL, NULL), ifelse(four, NULL, NULL))
spec_breaks <- NULL # ifelse(univariate, ifelse(sine, NULL, NULL), ifelse(four, NULL, NULL))
sens_plot <- plt_res(df_num$sens, "Sensitivity (%)", breaks = sens_breaks)
spec_plot <- plt_res(df_num$spec, "Specificity (%)", TeX("$\\textit{p}$"), breaks = spec_breaks)
time_plot <- plt_res(df_num$time, "Time (seconds)", log=T)

title <- "Results: $\\textit{q}="
plotname <-"q"
if (univariate){
  title <- paste0(title, "1$, ")
  plotname <- paste0(plotname, "1")
  if (sine){
    title <- paste0(title, "Non-Linear")
    plotname <- paste0(plotname, "NL")
  }else{
    title <- paste0(title, "Piece-Wise Linear")
    plotname <- paste0(plotname, "PWL")
  }
}else{
  if (four){
    title <- paste0(title, "4$")
    plotname <- paste0(plotname, "4")
  }else{
    title <- paste0(title, "2$")
    plotname <- paste0(plotname, "2")
  }
}
res_plots <- ggarrange(sens_plot, spec_plot, time_plot, ncol=3, common.legend = TRUE,  legend="bottom")
res_plots <- annotate_figure(res_plots, top = text_grob(TeX(title), size = 18, family = "Times"))
res_plots
ggsave(paste0("plots/", plotname, ".pdf"), height = 4, width = 14)

subgroups_list[[exper]] <- factor(Reduce(c, subgroups), levels = 2:7)

windowsFonts("Times" = windowsFont("Times"))
colors <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF") # https://nanx.me/ggsci/reference/pal_nejm.html
# plots[[exper]] <- list(x = factor(subgroups, levels = 2:6), plot = NULL)
plots <- list(ggplot() +
                geom_bar(aes(x = subgroups_list[["cont_cov_dep_"]]),
                         color = "black", fill = colors[1]) +
                theme_pubclean() +
                theme(text = element_text(family = "Times", size = 18),
                      plot.title = element_text(hjust = 0.5)) +
                ggtitle(TeX(paste0("Optimal $\\textit{K}, \\textit{q}=1$"))) +
                scale_y_continuous(breaks = seq(0, 180, 30), limits = c(0, 180)) +
                scale_x_discrete(drop=FALSE) +
                labs(x = TeX("$\\textit{K}$")),
              ggplot() +
                geom_bar(aes(x = subgroups_list[["cont_multi_cov_dep_"]]),
                         color = "black", fill = colors[2]) +
                theme_pubclean() +
                theme(text = element_text(family = "Times", size = 18),
                      plot.title = element_text(hjust = 0.5)) +
                ggtitle(TeX(paste0("Optimal $\\textit{K}, \\textit{q}=2$"))) +
                scale_y_continuous(breaks = seq(0, 180, 30), limits = c(0, 180)) +
                scale_x_discrete(drop=FALSE) +
                labs(x = TeX("$\\textit{K}$")),
              ggplot() +
                geom_bar(aes(x = subgroups_list[["cont_4_cov_dep_"]]),
                         color = "black", fill = colors[3]) +
                theme_pubclean() +
                theme(text = element_text(family = "Times", size = 18),
                      plot.title = element_text(hjust = 0.5)) +
                ggtitle(TeX(paste0("Optimal $\\textit{K}, \\textit{q}=4$"))) +
                scale_y_continuous(breaks = seq(0, 180, 30), limits = c(0, 180)) +
                scale_x_discrete(drop=FALSE) +
                labs(x = TeX("$\\textit{K}$")),
              ggplot() +
                geom_bar(aes(x = subgroups_list[["cont_cov_dep_sine_"]]),
                         color = "black", fill = colors[4]) +
                theme_pubclean() +
                theme(text = element_text(family = "Times", size = 18),
                      plot.title = element_text(hjust = 0.5)) +
                ggtitle(TeX(paste0("Optimal $\\textit{K}, \\textit{q}=1$ (non-linear)"))) +
                scale_y_continuous(breaks = seq(0, 180, 30), limits = c(0, 180)) +
                scale_x_discrete(drop=FALSE) +
                labs(x = TeX("$\\textit{K}$")))
fig <- ggarrange(plotlist = plots, nrow = 4)
ggsave("plots/subgroup_bar_newcol.pdf", height = 12, width = 11)

# ------------------------------------------------------------------------------
# HP comparisons
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each

# results config; switches from comparisons with max_iter_grid=max_iter=100 to
# max_iter_grid=10
max_iter_grid <- F

if (max_iter_grid){
  # comparisons with max_iter_grid = max_iter (set max_iter_grid = T)
  exper <- "cont_cov_dep_"
  n <- 75
  samp_str <- paste0("_n1_", n, "_n2_", n, "_n3_", n)
  path <- "./experiments/z1_hp"
  files <- list.files(path)
  files <- files[grepl("model_average", files)]
  files <- c(files, list.files(paste0(path, "_full")))
}else{
  # default comparisons
  exper <- "cont_cov_dep_"
  n <- 75
  samp_str <- paste0("_n1_", n, "_n2_", n, "_n3_", n)
  path <- "./experiments/z1_hp"
  files <- c(list.files(path), list.files(substr(path, 1, nchar(path) - 3)))
}
ntrials <- 50
trial_str <- paste0("ntrials", ntrials, "_")
dims <- c(10, 25, 50, 100)
prec <- 2

df <- list(grid_search = NA, hybrid = NA, model_average = NA)
df <- replicate(length(dims), df, simplify = F)
names(df) <- dims
df_num <- list(sens=data.frame(), spec=data.frame(), time=data.frame())

p <- "10"
hp_method <- names(df[[p]])[1]
for (p in as.character(dims)){

  best_sens <- best_spec <- 0
  best_time <- Inf
  sens_method <- spec_method <- time_method <- NA

  for (hp_method in names(df[[p]])){

    # format file name
    exp_name <- paste0(exper, trial_str, "p", p, samp_str)

    # load results

    file_name <- files[grepl(exp_name, files) & grepl(hp_method, files) & endsWith(files, ".Rda")]

    if (length(file_name) != 1) stop(paste(length(file_name), "files found"))
    file_path <- file.path(path, file_name)
    if (!file.exists(file_path)){
      if (max_iter_grid){
        file_path <- file.path(paste0(path, "_full"), file_name)
      }else{
        file_path <- file.path(substr(path, 1, nchar(path) - 3), file_name)
      }
    }
    load(file_path)
    print(paste0("n: ", results$sample_data[1], ", p: ", results$sample_data[2]))
    results <- results[setdiff(names(results), "sample_data")]

    results <- lapply(results, `[`, "covdepGE")

    # process sensitivity results
    sens <- sapply(results, sapply, `[[`, "sens") * 100
    sens_mean <- mean(sens)
    if (sens_mean > best_sens){
      best_sens <- sens_mean
      sens_method <- hp_method
    }
    sens_mean_str <- sprintf(paste0("%.", prec, "f"), sens_mean)
    sens_sd <- sd(sens)
    sens_se <- sens_sd / sqrt(ntrials)
    sens_se_str  <- sprintf(paste0("%.", prec, "f"), sens_se)
    sens_str <- paste0(sens_mean_str, " (\\new{", sens_se_str, "})")

    # process specificity results
    spec <- sapply(results, sapply, `[[`, "spec") * 100
    spec_mean <- mean(spec)
    if (spec_mean > best_spec){
      best_spec <- spec_mean
      spec_method <- hp_method
    }
    spec_mean_str <- sprintf(paste0("%.", prec, "f"), spec_mean)
    spec_sd <-  sd(spec)
    spec_se <- spec_sd / sqrt(ntrials)
    spec_se_str  <- sprintf(paste0("%.", prec, "f"), spec_se)
    spec_str <- paste0(spec_mean_str, " (\\new{", spec_se_str, "})")

    # process time results
    time <- sapply(results, sapply, `[[`, "time")
    time_mean <- mean(time)
    if (time_mean < best_time){
      best_time <- time_mean
      time_method <- hp_method
    }
    time_mean_str <- sprintf(paste0("%.", prec, "f"), time_mean)
    time_sd <- sd(time)
    time_se <- time_sd / sqrt(ntrials)
    time_se_str  <- sprintf(paste0("%.", prec, "f"), time_se)
    time_str <- paste0(time_mean_str, " (\\new{", time_se_str, "})")

    # combine summary strings
    perf_str <- cbind(sens_str, spec_str, time_str)
    perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
    row.names(perf_str) <- row.names(spec)

    # combine summary numerics
    df_num[["sens"]] <- rbind.data.frame(df_num[["sens"]],data.frame(p=p, Method=hp_method, mean=sens_mean, se=sens_se))
    df_num[["spec"]] <- rbind.data.frame(df_num[["spec"]],data.frame(p=p, Method=hp_method, mean=spec_mean, se=spec_se))
    df_num[["time"]] <- rbind.data.frame(df_num[["time"]],data.frame(p=p, Method=hp_method, mean=time_mean, se=time_se))


    # create storage
    df_exp <- data.frame(p = p, method = hp_method,
                         sens = NA, spec = NA, time = NA)

    df_exp[ , c("sens", "spec", "time")] <- perf_str

    df[[p]][[hp_method]] <- df_exp
    rm("results")
  }
  best_sens <- sprintf(paste0("%.", prec, "f"), best_sens * 100)
  best_spec <- sprintf(paste0("%.", prec, "f"), best_spec * 100)
  best_time <- sprintf(paste0("%.", prec, "f"), best_time)
  df[[p]][[sens_method]]$sens <- gsub(best_sens, paste0("\\\\mathbf{", best_sens, "}"), df[[p]][[sens_method]]$sens)
  df[[p]][[spec_method]]$spec <- gsub(best_spec, paste0("\\\\mathbf{", best_spec, "}"), df[[p]][[spec_method]]$spec)
  df[[p]][[time_method]]$time <- gsub(best_time, paste0("\\\\mathbf{", best_time, "}"), df[[p]][[time_method]]$time)
}

df <- lapply(df, function(x) Reduce(rbind, x))
df <- Reduce(rbind, df)
df$method <- gsub("_", "\\\\_", df$method)
df$method <- paste0("\\texttt{", df$method, "}")
colnames(df) <- c("$p$", "HP Method", "Sensitivity$(\\%)$", "Specificity$(\\%)$", "Time (seconds)")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

df_num

for (metric in names(df_num)){
  df_num[[metric]]$p <- factor(df_num[[metric]]$p, levels=c(10, 25, 50, 100))
  df_num[[metric]]$Method <- gsub("_sortZ", "_time", df_num[[metric]]$Method)
  # df_num[[metric]]$Method <- paste0("\\textit{", df_num[[metric]]$Method, "}")
  levels <- unique(df_num[[metric]]$Method)[order(unique(df_num[[metric]]$Method))]
  colors <- ggsci::pal_nejm()
  pal_colors <- colors(3)
  df_num[[metric]]$Method <- factor(df_num[[metric]]$Method, levels = levels)
  dodge <- unique(as.numeric(df_num[[metric]]$Method)/20)-mean(unique(as.numeric(df_num[[metric]]$Method)/20))
  df_num[[metric]]$p2 <- as.numeric(df_num[[metric]]$p) + dodge
}

plt_res <- function(data, ylab, xlab="", log=F, breaks=NULL, sz=0.75){
  plt <- ggplot(data, aes(x=p2, y=mean, group=Method, color=Method)) +
    theme_pubclean() + scale_color_manual(values = pal_colors, drop=FALSE) +
    theme(plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          text = element_text(family = "Times", size=18),
          legend.title=element_blank(), legend.key=element_blank(),
          legend.text=element_text(family="Courier", size=18)) +
    labs(x=xlab, y=ylab) +
    lapply(rev(levels), function(method) {
      list(
        geom_errorbar(data = ~ subset(., Method == method),
                      aes(ymin = mean - ifelse((log & mean - 2 * se < 0), se, 2 * se), ymax = mean + 2 * se),
                      width = 0.2, size = sz),
        geom_point(data = ~ subset(., Method == method), size = 2, aes()),
        geom_line(data = ~ subset(., Method == method), size = sz)
      )
    }) + scale_x_continuous(labels=c(10, 25, 50, 100), breaks = 1:4)

  if (log){
    plt <- plt + scale_y_continuous(trans='log10')
  }
  if (!is.null(breaks)){
    plt <- plt + scale_y_continuous(breaks = breaks)
  }
  plt
}

sens_plot <- plt_res(df_num$sens, "Sensitivity (%)")
spec_plot <- plt_res(df_num$spec, "Specificity (%)", TeX("$\\textit{p}$"))
time_plot <- plt_res(df_num$time, "Time (seconds)", log=T)

title <- "Results: Hyperparameter Specification Strategy Comparison"
plotname <-"hp_comp_plt"
res_plots <- ggarrange(sens_plot, spec_plot, time_plot, ncol=3, common.legend = TRUE,  legend="bottom")
res_plots <- annotate_figure(res_plots, top = text_grob(TeX(title), size = 18, family = "Times"))
res_plots
ggsave(paste0("plots/", plotname, ".pdf"), height = 4, width = 14)
