setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)
sf <- c(0.8, 0.99)
b <- qnorm(sf, 0, 1)
si <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)



file_names <- list.files("simulate_crosses_best_parents/offspring_family_info_true/")
info_df <- sapply(file_names, FUN=function(x){
  strsplit(x, "_")[[1]][1:5]
})
info_df <- t(info_df)
colnames(info_df) <- c("h2s", "effective_marker_sizes", "trait_number", "parent1", "parent2")
info_df <- as.data.frame(info_df)

info_df$h2s <- as.numeric(info_df$h2s)
info_df$effective_marker_sizes <- factor(info_df$effective_marker_sizes, 
                                levels=as.factor(effective_marker_sizes))
info_df$trait_number <- as.numeric(info_df$trait_number)
info_df$family <- paste(info_df$parent1, info_df$parent2, sep="_")
info_df <- info_df[, c("family", "h2s", "effective_marker_sizes", "trait_number")]

info_df$BV_mean <- NA
info_df$BV_sd <- NA

for (i in 1:nrow(info_df)){
  offspring = readRDS(paste("simulate_crosses_best_parents/offspring_family_info_true/", 
                            rownames(info_df)[i], sep=""))
  info_df$BV_mean[i] = offspring[[2]]
  info_df$BV_sd[i] = sqrt(offspring[[3]])
  
}

# file_names_pred <- list.files("simulate_crosses_best_parents/offspring_family_info_pred/")

info_df$predY_RR_mean <- NA
info_df$predY_RR_sd <- NA
info_df$predY_mean <- NA
info_df$predY_sd <- NA
for (i in 1:nrow(info_df)){
  offspring = readRDS(paste("simulate_crosses_best_parents/offspring_family_info_pred/", 
                            rownames(info_df)[i], sep=""))
  info_df[i, c("predY_RR_mean", "predY_RR_sd", "predY_mean", "predY_sd")] = 
    c(offspring[[2]], sqrt(offspring[[3]]), offspring[[5]], sqrt(offspring[[6]]))
}

saveRDS(info_df, "view_usefulness_best_parents/info_df.rds")
info_df <- readRDS("view_usefulness_best_parents/info_df.rds")



results_use <- info_df[rep(1:nrow(info_df), 2), ]
results_use$si <- rep(si, each=nrow(info_df))
results_use$BV_use <- results_use$BV_mean + results_use$si * results_use$BV_sd
results_use$predY_RR_use <- results_use$predY_RR_mean + results_use$si * results_use$predY_RR_sd
results_use$predY_use <- results_use$predY_mean + results_use$si * results_use$predY_sd

# 2*3*5*20
results_cor_use <- matrix(NA, nrow=600, ncol=9)
results_cor_use <- as.data.frame(results_cor_use)
colnames(results_cor_use) <- c("si", "h2s", "effective_marker_sizes", "trait_number", 
                               "use_mean_RR_cor", "use_use_RR_cor", "use_mean_cor", "use_use_cor", 
                               "realuse_realmean_cor")
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        # print((i-1)*5*20 + (j-1)*20 + k)
        A = results_use[results_use$si == si[h] & 
                          results_use$h2s == h2s[i] & 
                          results_use$effective_marker_sizes == effective_marker_sizes[j] & 
                          results_use$trait_number == k, ]
        results_cor_use[(h-1)*3*5*20 + (i-1)*5*20 + (j-1)*20 + k, ] = 
          c(si[h], h2s[i], effective_marker_sizes[j], k, 
            cor(A$BV_use, A$predY_RR_mean), cor(A$BV_use, A$predY_RR_use), 
            cor(A$BV_use, A$predY_mean), cor(A$BV_use, A$predY_use), 
            cor(A$BV_use, A$BV_mean))
      }
    }
  }
}
results_cor_use$effective_marker_sizes <- factor(results_cor_use$effective_marker_sizes, 
                                                 levels=as.factor(effective_marker_sizes))

saveRDS(results_cor_use, "view_usefulness_best_parents/results_cor_use.rds")
results_cor_use <- readRDS("view_usefulness_best_parents/results_cor_use.rds")

# 2*3*5
results_se_use <- matrix(NA, nrow=30, ncol=11)
results_se_use <- as.data.frame(results_se_use)
colnames(results_se_use) <- c("si", "h2s", "effective_marker_sizes", 
                              "use_mean_RR_mean", "use_mean_RR_se", "use_use_RR_mean", "use_use_RR_se", 
                              "use_mean_mean", "use_mean_se", "use_use_mean", "use_use_se")
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      A = results_cor_use[results_cor_use$si == si[h] & 
                            results_cor_use$h2s == h2s[i] & 
                            results_cor_use$effective_marker_sizes == effective_marker_sizes[j], ]
      results_se_use[(h-1)*3*5 + (i-1)*5 + j, ] = 
        c(si[h], h2s[i], effective_marker_sizes[j], 
          mean(A$use_mean_RR_cor), sd(A$use_mean_RR_cor)/sqrt(20), 
          mean(A$use_use_RR_cor), sd(A$use_use_RR_cor)/sqrt(20), 
          mean(A$use_mean_cor), sd(A$use_mean_cor)/sqrt(20), 
          mean(A$use_use_cor), sd(A$use_use_cor)/sqrt(20))
    }
  }
}
results_se_use$effective_marker_sizes <- factor(results_se_use$effective_marker_sizes, 
                                                levels=as.factor(effective_marker_sizes))

results_se_use$h2 <- paste("h^2 == ", results_se_use$h2s, sep="")
results_se_use$i <- paste("i == ", round(results_se_use$si, 2), sep="")
# pdf("view_usefulness/plots/test1.pdf")
p1 <- ggplot(results_se_use, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=use_mean_mean, colour="predicted from family mean"), size=0.75) + 
  geom_errorbar(aes(ymin=use_mean_mean-use_mean_se, ymax=use_mean_mean+use_mean_se, 
                    colour="predicted from family mean"), width=0.2) + 
  geom_point(aes(y=use_use_mean, colour="predicted from family usefulness"), size=0.75) + 
  geom_errorbar(aes(ymin=use_use_mean-use_use_se, ymax=use_use_mean+use_use_se, 
                    colour="predicted from family usefulness"), width=0.2) + 
  geom_line(aes(y=use_mean_mean, colour="predicted from family mean"), linewidth=0.5) + 
  geom_line(aes(y=use_use_mean, colour="predicted from family usefulness"), linewidth=0.5) + 
  facet_grid(i~h2, labeller = label_parsed) +
  # facet_wrap(~h2s, labeller = as_labeller(lbs, label_parsed)) + 
  xlab("number of causal loci") + 
  ylab("prediction accuracy\nof usefulness") + 
  ylim(-0.25, 0.75) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  guides(color=guide_legend(title="prediction method", nrow=2)) + 
  # ggtitle("accuracy of predicting family usefulness of BV, BayesC") + 
  theme_minimal_grid(font_size=7) +
  theme(legend.position="bottom") 
save_plot(paste("view_usefulness_best_parents/plots/", "fam_cor_best_parents.pdf", sep=""),
          plot_grid(p1),
          base_width=3.5, base_height=2.33)




















