# view_usefulness_best_parents_gametes.R



setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s                    <- c(0.8, 0.5, 0.2)
sf                     <- c(0.8, 0.99)
b                      <- qnorm(sf, 0, 1)
si                     <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)



info_df_best         <- readRDS("view_usefulness_best_parents/info_df.rds")
results_cor_use_best <- readRDS("view_usefulness_best_parents/results_cor_use.rds")



file_names <- list.files("simulate_crosses_best_parents_gametes/offspring_family_info_true/")
info_df    <- sapply(file_names, FUN=function(x){
  strsplit(x, "_")[[1]][1:5]
})
info_df           <- t(info_df)
colnames(info_df) <- c("h2s", "effective_marker_sizes", "trait_number", "parent1", "parent2")
info_df           <- as.data.frame(info_df)

info_df$h2s                    <- as.numeric(info_df$h2s)
info_df$effective_marker_sizes <- factor(info_df$effective_marker_sizes, 
                                         levels=as.factor(effective_marker_sizes))
info_df$trait_number           <- as.numeric(info_df$trait_number)
info_df$family                 <- paste(info_df$parent1, info_df$parent2, sep="_")

info_df$BV_mean <- NA
info_df$BV_sd <- NA
for (i in 1:nrow(info_df)){
  offspring = readRDS(paste("simulate_crosses_best_parents_gametes/offspring_family_info_true/", 
                            rownames(info_df)[i], sep=""))
  info_df$BV_mean[i] = offspring[[2]]
  info_df$BV_sd[i]   = sqrt(offspring[[3]])
}

info_df$predY_RR_mean <- NA
info_df$predY_RR_sd   <- NA
info_df$predY_mean    <- NA
info_df$predY_sd      <- NA
for (i in 1:nrow(info_df)){
  offspring = readRDS(paste("simulate_crosses_best_parents_gametes/offspring_family_info_pred/", 
                            rownames(info_df)[i], sep=""))
  info_df[i, c("predY_RR_mean", "predY_RR_sd", "predY_mean", "predY_sd")] = 
    c(offspring[[2]], sqrt(offspring[[3]]), offspring[[5]], sqrt(offspring[[6]]))
}
info_df$si <- NA
info_df    <- info_df[, c(13, 1:12)]
dim(info_df)
# [1] 3000   13
info_df    <- info_df[rep(1:nrow(info_df), 2), ]
info_df$si <- rep(round(si, 2), each=nrow(info_df)/2)
dim(info_df)
# [1] 6000   13

info_df$BV_use       <- info_df$BV_mean + info_df$si * info_df$BV_sd
info_df$predY_RR_use <- info_df$predY_RR_mean + info_df$si * info_df$predY_RR_sd
info_df$predY_use    <- info_df$predY_mean + info_df$si * info_df$predY_sd

saveRDS(info_df, "view_usefulness_best_parents_gametes/info_df.rds")
info_df <- readRDS("view_usefulness_best_parents_gametes/info_df.rds")





file_names <- list.files("simulate_crosses_best_parents_gametes_use/offspring_family_info_true/")
info_df_use <- sapply(file_names, FUN=function(x){
  strsplit(x, "_")[[1]][1:6]
})
info_df_use <- t(info_df_use)
colnames(info_df_use) <- c("si", "h2s", "effective_marker_sizes", "trait_number", "parent1", "parent2")
info_df_use <- as.data.frame(info_df_use)

info_df_use$si                     <- as.numeric(info_df_use$si)
info_df_use$h2s                    <- as.numeric(info_df_use$h2s)
info_df_use$effective_marker_sizes <- factor(info_df_use$effective_marker_sizes, 
                                         levels=as.factor(effective_marker_sizes))
info_df_use$trait_number           <- as.numeric(info_df_use$trait_number)
info_df_use$family                 <- paste(info_df_use$parent1, info_df_use$parent2, sep="_")

info_df_use$BV_mean <- NA
info_df_use$BV_sd   <- NA
for (i in 1:nrow(info_df_use)){
  offspring = readRDS(paste("simulate_crosses_best_parents_gametes_use/offspring_family_info_true/", 
                            rownames(info_df_use)[i], sep=""))
  info_df_use$BV_mean[i] = offspring[[2]]
  info_df_use$BV_sd[i]   = sqrt(offspring[[3]])
}

info_df_use$predY_RR_mean <- NA
info_df_use$predY_RR_sd   <- NA
info_df_use$predY_mean    <- NA
info_df_use$predY_sd      <- NA
for (i in 1:nrow(info_df_use)){
  offspring = readRDS(paste("simulate_crosses_best_parents_gametes_use/offspring_family_info_pred/", 
                            rownames(info_df_use)[i], sep=""))
  info_df_use[i, c("predY_RR_mean", "predY_RR_sd", "predY_mean", "predY_sd")] = 
    c(offspring[[2]], sqrt(offspring[[3]]), offspring[[5]], sqrt(offspring[[6]]))
}
dim(info_df_use)
# [1] 6000   13

info_df_use$BV_use       <- info_df_use$BV_mean + info_df_use$si * info_df_use$BV_sd
info_df_use$predY_RR_use <- info_df_use$predY_RR_mean + info_df_use$si * info_df_use$predY_RR_sd
info_df_use$predY_use    <- info_df_use$predY_mean + info_df_use$si * info_df_use$predY_sd

saveRDS(info_df_use, "view_usefulness_best_parents_gametes/info_df_use.rds")
info_df_use <- readRDS("view_usefulness_best_parents_gametes/info_df_use.rds")



# 2*3*5*20
results           <- matrix(NA, nrow=600, ncol=9)
results           <- as.data.frame(results)
colnames(results) <- c("si", "h2s", "effective_marker_sizes", "trait_number", 
                       "use_mean_RR_cor", "use_use_RR_cor", "use_mean_cor", "use_use_cor", 
                       "realuse_realmean_cor")
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        # print((i-1)*5*20 + (j-1)*20 + k)
        A = info_df[info_df$si == round(si[h], 2) & 
                      info_df$h2s == h2s[i] & 
                      info_df$effective_marker_sizes == effective_marker_sizes[j] & 
                      info_df$trait_number == k, ]
        results[(h-1)*3*5*20 + (i-1)*5*20 + (j-1)*20 + k, ] = 
          c(round(si[h], 2), h2s[i], effective_marker_sizes[j], k, 
            cor(A$BV_use, A$predY_RR_mean), cor(A$BV_use, A$predY_RR_use), 
            cor(A$BV_use, A$predY_mean), cor(A$BV_use, A$predY_use), 
            cor(A$BV_use, A$BV_mean))
      }
    }
  }
}
results$effective_marker_sizes <- factor(results$effective_marker_sizes, 
                                         levels=as.factor(effective_marker_sizes))

saveRDS(results, "view_usefulness_best_parents_gametes/results.rds")
results <- readRDS("view_usefulness_best_parents_gametes/results.rds")

results_use           <- matrix(NA, nrow=600, ncol=9)
results_use           <- as.data.frame(results_use)
colnames(results_use) <- c("si", "h2s", "effective_marker_sizes", "trait_number", 
                       "use_mean_RR_cor", "use_use_RR_cor", "use_mean_cor", "use_use_cor", 
                       "realuse_realmean_cor")
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        # print((i-1)*5*20 + (j-1)*20 + k)
        A = info_df_use[info_df_use$si == round(si[h], 2) & 
                          info_df_use$h2s == h2s[i] & 
                          info_df_use$effective_marker_sizes == effective_marker_sizes[j] & 
                          info_df_use$trait_number == k, ]
        results_use[(h-1)*3*5*20 + (i-1)*5*20 + (j-1)*20 + k, ] = 
          c(round(si[h], 2), h2s[i], effective_marker_sizes[j], k, 
            cor(A$BV_use, A$predY_RR_mean), cor(A$BV_use, A$predY_RR_use), 
            cor(A$BV_use, A$predY_mean), cor(A$BV_use, A$predY_use), 
            cor(A$BV_use, A$BV_mean))
      }
    }
  }
}
results_use$effective_marker_sizes <- factor(results_use$effective_marker_sizes, 
                                             levels=as.factor(effective_marker_sizes))

saveRDS(results_use, "view_usefulness_best_parents_gametes/results_use.rds")
results_use <- readRDS("view_usefulness_best_parents_gametes/results_use.rds")



results_se           <- matrix(NA, nrow=30, ncol=11)
results_se           <- as.data.frame(results_se)
colnames(results_se) <- c("si", "h2s", "effective_marker_sizes", 
                          "use_mean_RR_mean", "use_mean_RR_se", "use_use_RR_mean", "use_use_RR_se", 
                          "use_mean_mean", "use_mean_se", "use_use_mean", "use_use_se")
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      A = results[results$si == round(si[h], 2) & 
                            results$h2s == h2s[i] & 
                            results$effective_marker_sizes == effective_marker_sizes[j], ]
      results_se[(h-1)*3*5 + (i-1)*5 + j, ] = 
        c(round(si[h], 2), h2s[i], effective_marker_sizes[j], 
          mean(A$use_mean_RR_cor), sd(A$use_mean_RR_cor)/sqrt(20), 
          mean(A$use_use_RR_cor), sd(A$use_use_RR_cor)/sqrt(20), 
          mean(A$use_mean_cor), sd(A$use_mean_cor)/sqrt(20), 
          mean(A$use_use_cor), sd(A$use_use_cor)/sqrt(20))
    }
  }
}
results_se$effective_marker_sizes <- factor(results_se$effective_marker_sizes, 
                                            levels=as.factor(effective_marker_sizes))

results_se_use           <- matrix(NA, nrow=30, ncol=11)
results_se_use           <- as.data.frame(results_se_use)
colnames(results_se_use) <- c("si", "h2s", "effective_marker_sizes", 
                          "use_mean_RR_mean", "use_mean_RR_se", "use_use_RR_mean", "use_use_RR_se", 
                          "use_mean_mean", "use_mean_se", "use_use_mean", "use_use_se")
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      A = results_use[results_use$si == round(si[h], 2) & 
                        results_use$h2s == h2s[i] & 
                        results_use$effective_marker_sizes == effective_marker_sizes[j], ]
      results_se_use[(h-1)*3*5 + (i-1)*5 + j, ] = 
        c(round(si[h], 2), h2s[i], effective_marker_sizes[j], 
          mean(A$use_mean_RR_cor), sd(A$use_mean_RR_cor)/sqrt(20), 
          mean(A$use_use_RR_cor), sd(A$use_use_RR_cor)/sqrt(20), 
          mean(A$use_mean_cor), sd(A$use_mean_cor)/sqrt(20), 
          mean(A$use_use_cor), sd(A$use_use_cor)/sqrt(20))
    }
  }
}
results_se_use$effective_marker_sizes <- factor(results_se_use$effective_marker_sizes, 
                                            levels=as.factor(effective_marker_sizes))


# this is only on the BayesC model
results_long                 <- results_se[rep(1:nrow(results_se), 2), 1:3]
results_long$prediction_type <- rep(c("predicted from family mean", "predicted from family usefulness"), 
                                    each=nrow(results_se))
results_long$parents_type    <- "selected based on predicted mean"
results_long$mean            <- c(results_se$use_mean_mean, results_se$use_use_mean)
results_long$se              <- c(results_se$use_mean_se, results_se$use_mean_se)

results_long_use                 <- results_se_use[rep(1:nrow(results_se_use), 2), 1:3]
results_long_use$prediction_type <- rep(c("predicted from family mean", "predicted from family usefulness"), 
                                    each=nrow(results_se_use))
results_long_use$parents_type    <- "selected based on predicted usefulness"
results_long_use$mean            <- c(results_se_use$use_mean_mean, results_se_use$use_use_mean)
results_long_use$se              <- c(results_se_use$use_mean_se, results_se_use$use_mean_se)

results_long      <- rbind(results_long, results_long_use)
results_long$type <- paste(results_long$prediction_type, results_long$parents_type, sep=", ")
results_long$i    <- paste("i == ", results_long$si, sep="")
results_long$h2   <- paste("h^2 == ", results_long$h2s, sep="")



p1 <- ggplot(results_long, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=mean, color=type)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color=type), width=0.2) + 
  geom_line(aes(y=mean, color=type), linewidth=0.5) + 
  facet_grid(i~h2, labeller = label_parsed) +
  # facet_wrap(~h2s, labeller = as_labeller(lbs, label_parsed)) + 
  xlab("number of causal loci") + 
  ylab("prediction accuracy\nof usefulness") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  guides(color=guide_legend(title="prediction method\nand parents type", ncol=1)) + 
  # ggtitle("accuracy of predicting family usefulness of BV, BayesC") + 
  theme_minimal_grid(font_size=10) +
  theme(legend.position="bottom") 
save_plot(paste("view_usefulness_best_parents_gametes/plots/", "fam_cor_best_parents_gametes.pdf", sep=""),
          plot_grid(p1),
          base_width=6.5, base_height=5)



