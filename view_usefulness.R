setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)
sf <- c(0.8, 0.99)
b <- qnorm(sf, 0, 1)
si <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)

# file_names_BV <- list.files("simulate_phenotypes_crosses/", pattern="BV.rds")
# file_names_BV_fammean <- list.files("simulate_phenotypes_crosses/", pattern="BV_fammean.rds")
# file_names_BV_famvar <- list.files("simulate_phenotypes_crosses/", pattern="BV_famvar.rds")
# file_names_predY_RR_fammean <- list.files("simulate_phenotypes_crosses/", 
#                                           pattern="predY_RR_fammean.rds")
# file_names_predY_RR_famvar <- list.files("simulate_phenotypes_crosses/", 
#                                          pattern="predY_RR_famvar.rds")
# file_names_predY_fammean <- list.files("simulate_phenotypes_crosses/", 
#                                        pattern="predY_fammean.rds")
# file_names_predY_famvar <- list.files("simulate_phenotypes_crosses/", 
#                                       pattern="predY_famvar.rds")
file_names_results <- list.files("simulate_phenotypes_crosses/", pattern="_result_df.rds")
file_names <- gsub("_result_df.rds", "", file_names_results)
parents_names <- matrix(NA, nrow=length(file_names), ncol=2)
parents_names[, 1] <- gsub("_.+", "", file_names)
parents_names[, 2] <- gsub(".+_", "", file_names)
rownames(parents_names) <- file_names

# # make sure parents_names is in the right order
# names_parents <-  read.delim("simulate_crosses/file_names.txt", header = F)
# names_parents <- gsub(".rds", "", gsub("simulate_crosses/offspring_list_", "", names_parents[, 1]))
# all(names_parents == rownames(parents_names))
# # [1] TRUE

Z <- readRDS("simulate_phenotypes/Z.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")

parents_Z_noneffective <- vector(mode="list", length=length(effective_marker_indices))
names(parents_Z_noneffective) <- effective_marker_indices
parents_Z_noneffective <- lapply(parents_Z_noneffective, FUN=function(x){
  vector(mode="list", length=20)})
for (j in 1:length(parents_Z_noneffective)){
  parents_Z_noneffective[[j]] = lapply(parents_Z_noneffective[[j]], FUN=function(x){
    vector(mode="list", length=500)})
  for (k in 1:20){
    for (h in 1:500){
      parents_Z_noneffective[[j]][[k]][[h]] = 
        Z[parents_names[(k-1)*500 + h, ], -effective_marker_indices[[j]][, k]]
    }
  }
}
parents_heter_sites_noncausal <- vector(mode="list", length=length(effective_marker_indices))
names(parents_heter_sites_noncausal) <- effective_marker_sizes
for (j in 1:length(parents_heter_sites_noncausal)){
  parents_heter_sites_noncausal[[j]] = sapply(parents_Z_noneffective[[j]][[1]], FUN=function(x){
    sum(x==0)/length(x)
  })
  for (k in 2:20){
    parents_heter_sites_noncausal[[j]] = cbind(parents_heter_sites_noncausal[[j]], 
                                               sapply(parents_Z_noneffective[[j]][[k]], FUN=function(x){
                                                 sum(x==0)/length(x)
                                               }))
  }
}





# BV <- lapply(file_names_BV, FUN=function(x){
#   readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
# })
# for (i in 1:length(BV)){
#   BV[[i]] = BV[[i]][[1]]
# }

temp <- lapply(file_names_results, FUN=function(x) {
  readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
})
for (i in 1:length(temp)){
  temp[[i]]$BV_sd = sqrt(temp[[i]]$BV_var)
  temp[[i]]$predY_RR_sd = sqrt(temp[[i]]$predY_RR_var)
  temp[[i]]$predY_sd = sqrt(temp[[i]]$predY_var)
  temp[[i]] = temp[[i]][, -c(6, 8, 10)]
}



# organize the data
# in "results", each row is a family
results = matrix(NA, nrow=15*10000, ncol=10)
results <- as.data.frame(results)
for (i in 1:length(temp)){
  # if(i %% 500 == 0){
  #   print(i)
  # }
  results[((i-1)*15 + 1) : ((i-1)*15 + 15), ] = temp[[i]]
}
colnames(results) <- colnames(temp[[1]])
results$effective_marker_sizes <- factor(results$effective_marker_sizes, labels = effective_marker_sizes)
saveRDS(results, "view_usefulness/results.rds")
results <- readRDS("view_usefulness/results.rds")



# 300 = 5*3*20
results_cor <- matrix(NA, nrow=300, ncol=8)
results_cor <- as.data.frame(results_cor)
colnames(results_cor) <- c("h2s", "effective_marker_sizes", "trait_number", 
                           "mean_RR_cor", "mean_cor", "sd_RR_cor", "sd_cor", "sd_het_cor")
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      # print((i-1)*5*20 + (j-1)*20 + k)
      A = results[results$h2s == h2s[i] & 
                    results$effective_marker_sizes == effective_marker_sizes[j] & 
                    results$trait_number == k, ]
      results_cor[(i-1)*5*20 + (j-1)*20 + k, ] = 
        c(h2s[i], effective_marker_sizes[j], k, 
          cor(A$BV_mean, A$predY_RR_mean), cor(A$BV_mean, A$predY_mean), 
          cor(A$BV_sd, A$predY_RR_sd), cor(A$BV_sd, A$predY_sd), 
          cor(A$BV_sd, parents_heter_sites_noncausal[[j]][, k]))
    }
  }
}
results_cor$effective_marker_sizes <- factor(results_cor$effective_marker_sizes, 
                                             levels=as.factor(effective_marker_sizes))



# 15 = 3*5
results_se <- matrix(NA, nrow=15, ncol=12)
results_se <- as.data.frame(results_se)
colnames(results_se) <- c("h2s", "effective_marker_sizes", 
                               "mean_RR_mean", "mean_RR_se", "mean_mean", "mean_se", 
                               "sd_RR_mean", "sd_RR_se", "sd_mean", "sd_se", 
                               "sd_het_mean", "sd_het_se")
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    A = results_cor[results_cor$h2s == h2s[i] & 
                      results_cor$effective_marker_sizes == effective_marker_sizes[j], ]
    results_se[(i-1)*5 + j, ] = 
      c(h2s[i], effective_marker_sizes[j], 
        mean(A$mean_RR_cor), sd(A$mean_RR_cor)/sqrt(20), mean(A$mean_cor), sd(A$mean_cor)/sqrt(20), 
        mean(A$sd_RR_cor), sd(A$sd_RR_cor)/sqrt(20), mean(A$sd_cor), sd(A$sd_cor)/sqrt(20), 
        mean(A$sd_het_cor), sd(A$sd_het_cor)/sqrt(20))
  }
}
results_se$effective_marker_sizes <- factor(results_se$effective_marker_sizes, 
                                             levels=as.factor(effective_marker_sizes))



# plotting the prediction accuracy of family mean and variance 
results_se$h2 <- paste("h^2 == ", results_se$h2s, sep="")
# pdf("view_usefulness/plots/test1.pdf")
p1 <- ggplot(results_se, aes(as.numeric(effective_marker_sizes), mean_mean)) + 
  geom_point(color="blue", size=0.75) + 
  geom_errorbar(aes(ymin=mean_mean-mean_se, ymax=mean_mean+mean_se), width=0.2, color="blue") + 
  geom_line(color="blue", linewidth=0.5) +
  facet_wrap(~h2, labeller = label_parsed) + 
  xlab("number of causal loci") + 
  ylab("prediction accuracy\nof mean") + 
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  # ggtitle("accuracy of predicting family mean of BV, BayesC") + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# dev.off()

# save_plot(paste("view_usefulness/plots/", "fammean_cor_2.pdf", sep=""),
#           p1, 
#           base_width=6.5, base_height=2.17)

# pdf("view_usefulness/plots/test1.pdf")
p2 <- ggplot(results_se, aes(as.numeric(effective_marker_sizes))) +
  geom_point(aes(y=sd_mean, colour="prediction from model"), size=0.75) +
  geom_errorbar(aes(ymin=sd_mean-sd_se, ymax=sd_mean+sd_se, 
                    colour="prediction from model"), width=0.2) + 
  geom_line(aes(as.numeric(effective_marker_sizes), y=sd_mean, colour="prediction from model"), 
              linewidth=0.5) +
  geom_point(aes(y=sd_het_mean, colour="prediction from parental heterozygosity"), size=0.75) +
  geom_errorbar(aes(ymin=sd_het_mean-sd_het_se, ymax=sd_het_mean+sd_het_se, 
                    colour="prediction from parental heterozygosity"), width=0.2) + 
  geom_line(aes(as.numeric(effective_marker_sizes), y=sd_het_mean, 
                colour="prediction from parental heterozygosity"), 
              linewidth=0.5) + 
  facet_wrap(~h2, labeller = label_parsed) + 
  # ylim(-0.25, 1) +
  xlab("number of causal loci") +
  ylab("prediction accunracy\nof standard deviation") +
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(strip.text.x = element_blank(), 
        legend.position = "bottom"
        , axis.title.y = element_text(hjust = 1.5)
        ) + 
  guides(colour=guide_legend(title="prediction method", nrow=2))
#       axis.title.x=element_blank(),
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())
# theme(strip.background = element_blank())
# dev.off()

# prow <- plot_grid(p1, 
#                   p2, 
#                   labels="auto", ncol=1, rel_heights = c(0.85, 1))

save_plot(paste("view_usefulness/plots/", "fammean_famsd_cor.pdf", sep=""), 
          plot_grid(p1, p2, labels="auto", label_size = 7, ncol=1, rel_heights = c(0.725, 1)), 
          base_width=3.5, base_height = 2.33)

p8 <- ggplot(results_se, aes(as.numeric(effective_marker_sizes), mean_RR_mean)) + 
  geom_point(color="blue", size=0.75) + 
  geom_errorbar(aes(ymin=mean_RR_mean-mean_RR_se, ymax=mean_RR_mean+mean_RR_se), width=0.2, color="blue") + 
  geom_line(color="blue", linewidth=0.5) +
  facet_wrap(~h2, labeller = label_parsed) + 
  xlab("number of causal loci") + 
  ylab("prediction accuracy\nof mean") + 
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  # ggtitle("accuracy of predicting family mean of BV, BayesC") + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p9 <- ggplot(results_se, aes(as.numeric(effective_marker_sizes))) +
  geom_point(aes(y=sd_RR_mean, colour="prediction from model"), size=0.75) +
  geom_errorbar(aes(ymin=sd_RR_mean-sd_RR_se, ymax=sd_RR_mean+sd_RR_se, 
                    colour="prediction from model"), width=0.2) + 
  geom_line(aes(as.numeric(effective_marker_sizes), y=sd_RR_mean, colour="prediction from model"), 
            linewidth=0.5) +
  geom_point(aes(y=sd_het_mean, colour="prediction from parental heterozygosity"), size=0.75) +
  geom_errorbar(aes(ymin=sd_het_mean-sd_het_se, ymax=sd_het_mean+sd_het_se, 
                    colour="prediction from parental heterozygosity"), width=0.2) + 
  geom_line(aes(as.numeric(effective_marker_sizes), y=sd_het_mean, 
                colour="prediction from parental heterozygosity"), 
            linewidth=0.5) + 
  facet_wrap(~h2, labeller = label_parsed) + 
  # ylim(-0.25, 1) +
  xlab("number of causal loci") +
  ylab("prediction accuracy\nof standard deviation") +
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(strip.text.x = element_blank(), 
        legend.position = "bottom"
        , axis.title.y = element_text(hjust = 1.5)
        ) + 
  guides(colour=guide_legend(title="prediction method", nrow=2))
save_plot(paste("view_usefulness/plots/", "fammean_famsd_cor_RR.pdf", sep=""), 
          plot_grid(p8, p9, labels="auto", label_size = 7, ncol=1, rel_heights = c(0.725, 1)), 
          base_width=3.5, base_height = 2.33)



# create a results table with usefulness
# look at the correlation between true usefulness and predicted usefulness, 
# and between true usefulness and predicted mean
results_use <- results[rep(1:nrow(results), 2), ]
results_use$si <- rep(si, each=nrow(results))
results_use <- results_use[, c(1, 11, 2:10)]
results_use$BV_use <- results_use$BV_mean + results_use$si * results_use$BV_sd
results_use$predY_RR_use <- results_use$predY_RR_mean + results_use$si * results_use$predY_RR_sd
results_use$predY_use <- results_use$predY_mean + results_use$si * results_use$predY_sd

# 2*5*3*20
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

# saveRDS(results_cor_use, "view_usefulness/results_cor_use.rds")
# results_cor_use <- readRDS("view_usefulness/results_cor_use.rds")

# 30 = 2*3*5
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
p3 <- ggplot(results_se_use, aes(as.numeric(effective_marker_sizes))) + 
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
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  guides(color=guide_legend(title="prediction method", nrow=2)) + 
  # ggtitle("accuracy of predicting family usefulness of BV, BayesC") + 
  theme_minimal_grid(font_size=7) +
  theme(legend.position="bottom") 
save_plot(paste("view_usefulness/plots/", "famuse_cor.pdf", sep=""),
          plot_grid(p3),
          base_width=3.5, base_height=2.33)
# dev.off()

p10 <- ggplot(results_se_use, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=use_mean_RR_mean, colour="predicted from family mean"), size=0.75) + 
  geom_errorbar(aes(ymin=use_mean_RR_mean-use_mean_RR_se, ymax=use_mean_RR_mean+use_mean_RR_se, 
                    colour="predicted from family mean"), width=0.2) + 
  geom_point(aes(y=use_use_RR_mean, colour="predicted from family usefulness"), size=0.75) + 
  geom_errorbar(aes(ymin=use_use_RR_mean-use_use_RR_se, ymax=use_use_RR_mean+use_use_RR_se, 
                    colour="predicted from family usefulness"), width=0.2) + 
  geom_line(aes(y=use_mean_RR_mean, colour="predicted from family mean"), linewidth=0.5) + 
  geom_line(aes(y=use_use_RR_mean, colour="predicted from family usefulness"), linewidth=0.5) + 
  facet_grid(i~h2, labeller = label_parsed) +
  # facet_wrap(~h2s, labeller = as_labeller(lbs, label_parsed)) + 
  xlab("number of causal loci") + 
  ylab("prediction accuracy\nof usefulness") + 
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  guides(color=guide_legend(title="prediction method", nrow=2)) + 
  # ggtitle("accuracy of predicting family usefulness of BV, BayesC") + 
  theme_minimal_grid(font_size=7) +
  theme(legend.position="bottom") 
save_plot(paste("view_usefulness/plots/", "famuse_cor_RR.pdf", sep=""),
          plot_grid(p10),
          base_width=3.5, base_height=2.33)



# look at the correlation between true mean and true usefulness
# look at the variance of mean and variance of i*sd
# 2*5*20
results_cor_usemean_varmean_varsisd <- matrix(NA, nrow=200, ncol=7)
results_cor_usemean_varmean_varsisd <- as.data.frame(results_cor_usemean_varmean_varsisd)
colnames(results_cor_usemean_varmean_varsisd) <- c("si", "effective_marker_sizes", "trait_number", 
                               "use_mean_cor", "var_mean", "var_sd", "var_si_sd")
for (h in 1:length(si)){
  for (i in 1){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        # print((i-1)*5*20 + (j-1)*20 + k)
        A = results_use[results_use$si == si[h] & 
                          results_use$h2s == h2s[i] & 
                          results_use$effective_marker_sizes == effective_marker_sizes[j] & 
                          results_use$trait_number == k, ]
        results_cor_usemean_varmean_varsisd[(h-1)*5*20 + (j-1)*20 + k, ] = 
          c(si[h], effective_marker_sizes[j], k, 
            cor(A$BV_use, A$BV_mean), var(A$BV_mean), var(A$BV_sd), var(A$BV_sd * A$si))
      }
    }
  }
}

results_cor_usemean_varmean_varsisd$t_simulation <-
  results_cor_usemean_varmean_varsisd$var_sd / results_cor_usemean_varmean_varsisd$var_mean
results_cor_usemean_varmean_varsisd$t_analytical <-
  1 / (4*as.numeric(as.character(results_cor_usemean_varmean_varsisd$effective_marker_sizes)))
# results_cor_usemean_varmean_varsisd$L <- 
#   1 / (4*results_cor_usemean_varmean_varsisd$t_simulation)

results_cor_usemean_varmean_varsisd$effective_marker_sizes <- 
  factor(results_cor_usemean_varmean_varsisd$effective_marker_sizes, 
         levels=as.factor(effective_marker_sizes))
results_cor_usemean_varmean_varsisd$i <- 
  paste("i == ", round(results_cor_usemean_varmean_varsisd$si, 2), sep="")

# 2*5
results_se_usemean <- matrix(NA, nrow=10, ncol=4)
results_se_usemean <- as.data.frame(results_se_usemean)
colnames(results_se_usemean) <- c("si", "effective_marker_sizes", 
                              "use_mean_mean", "use_mean_se")
for (h in 1:length(si)){
  for (i in 1){
    for (j in 1:length(effective_marker_sizes)){
      A = results_cor_use[results_cor_use$si == si[h] & 
                            results_cor_use$h2s == h2s[i] & 
                            results_cor_use$effective_marker_sizes == effective_marker_sizes[j], ]
      results_se_usemean[(h-1)*5 + j, ] = 
        c(si[h], effective_marker_sizes[j], 
          mean(A$use_mean_cor), sd(A$use_mean_cor)/sqrt(20))
    }
  }
}
results_se_usemean$effective_marker_sizes <- factor(results_se_usemean$effective_marker_sizes, 
                                                levels=as.factor(effective_marker_sizes))

results_var <- matrix(NA, nrow=5, ncol=4)
results_var <- as.data.frame(results_var)
colnames(results_var) <- c("effective_marker_sizes", "t_type", "t_mean", "t_se")

for (j in 1:length(effective_marker_sizes)){
  for (h in 1:2){
    A = results_cor_usemean_varmean_varsisd[
      results_cor_usemean_varmean_varsisd$si == si[1] & 
        results_cor_usemean_varmean_varsisd$effective_marker_sizes == effective_marker_sizes[j], ]
    if (h == 1){
      results_var[(h-1)*5 + j, 1:4] = 
        c(effective_marker_sizes[j], "simulated offspring", 
          mean(A$t_simulation), sd(A$t_simulation)/sqrt(20))
    } else {
      results_var[(h-1)*5 + j, 1:4] = 
        c(effective_marker_sizes[j], "the number of causal loci", 
          mean(A$t_analytical), sd(A$t_analytical)/sqrt(20))
    }
  }
}

# for (h in 1:2){
#   for (i in 1){
#     for (j in 1:length(effective_marker_sizes)){
#       A = results_cor_usemean_varmean_varsisd[
#         results_cor_usemean_varmean_varsisd$si == si[1] & 
#           results_cor_usemean_varmean_varsisd$effective_marker_sizes == effective_marker_sizes[j], ]
#       results_var[j, 1:4] = 
#         c(effective_marker_sizes[j], 
#           mean(A$t_simulation), sd(A$t_simulation)/sqrt(20), 
#           mean(A$t_analytical))
#     }
#   }
# }
results_var$effective_marker_sizes <- 
  factor(results_var$effective_marker_sizes, 
         levels=as.factor(effective_marker_sizes))
results_var$t_mean <- as.numeric(results_var$t_mean)
results_var$t_se <- as.numeric(results_var$t_se)



results_se_usemean$i <- as.character(round(results_se_usemean$si, 2))
# pdf("view_usefulness/plots/test1.pdf")
p4 <-
  ggplot(results_se_usemean, aes(x=as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=use_mean_mean, color=i), size=0.75) + 
  geom_errorbar(aes(ymin=use_mean_mean-use_mean_se, ymax=use_mean_mean+use_mean_se, color=i), width=0.2) + 
  geom_line(aes(y=use_mean_mean, color=i), linewidth=0.5) + 
  xlab("number of causal loci") + 
  ylab("correlation between BV\nmean and usefulness") + 
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) +
  scale_colour_manual(values=c("blue", "gold3")) + 
  guides(color=guide_legend(title="selection\nintensity")) + 
  # scale_color_discrete(name="selection intensity", labels = parse(text=unique(results_se_usemean$i))) + 
  theme_minimal_grid(font_size=7) + 
  theme(
    # legend.position="top",
        # axis.title.y = element_text(hjust = 0.2),
        # axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank()
    )
# dev.off()
pdf("view_usefulness/plots/t_value.pdf", width=3.5, height=3.5)
# p5 <-
  ggplot(results_var, aes(as.numeric(effective_marker_sizes)
                          # , color=t_type, group=interaction(effective_marker_sizes, t_type)
                          )) + 
  geom_point(aes(y=t_mean, color=t_type), size=0.75) + 
  geom_errorbar(aes(ymin=t_mean-t_se, ymax=t_mean+t_se, color=t_type), width=0.2) + 
  geom_line(aes(y=t_mean, color=t_type), linewidth=0.5) + 
  # geom_boxplot(outlier.size = 0.75, size=0.4) + 
  # facet_wrap(~i, labeller = label_parsed) + 
  # scale_y_log10() + 
  xlab("number of causal loci") + 
  ylab("t") + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  # scale_y_continuous(trans="log10") +
  # scale_colour_manual(values=c("gold3", "blue")) + 
  guides(color=guide_legend(title="t value calculated\nbased on", nrow=2)) +
  theme_minimal_grid(font_size=7) + 
  theme(legend.position="bottom")
# p5 <- 
#   ggplot(results_var, aes(effective_marker_sizes, variance, color=variance_type)) + 
#   geom_boxplot(position = "dodge") + 
#   facet_wrap(~i, labeller = label_parsed) + 
#   scale_y_log10() + 
#   xlab("number of causal loci") + 
#   ylab("variance of mean or i*sd") + 
#   scale_colour_manual(values=c("gold3", "blue")) + 
#   guides(color=guide_legend(title="variance type")) + 
#   theme_minimal_grid(font_size=7) 
dev.off()
# pdf("view_usefulness/plots/test1.pdf")
# p5 <-
#   ggplot(results_cor_usemean_varmean_varsisd[
#   results_cor_usemean_varmean_varsisd$effective_marker_sizes=="4" & 
#     results_cor_usemean_varmean_varsisd$i == "i == 2.67", ], 
#              aes(var_mean, var_si_sd)) + 
#   geom_point() + 
#   xlim(0, 2.5) + 
#   ylim(0, 2.5) + 
#   geom_abline(slope=1, intercept=0) + 
#   # facet_wrap(~i, labeller = label_parsed) + 
#   xlab("variance of mean") + 
#   ylab("variance of i*sd") +
#   theme_minimal_grid(font_size=8) + 
#   theme(plot.title=element_text(size=8)) +
#   ggtitle("number of causal loci: 4")
# # dev.off()
# # pdf("view_usefulness/plots/test1.pdf")
# p6 <-
# ggplot(results_cor_usemean_varmean_varsisd[
#   results_cor_usemean_varmean_varsisd$effective_marker_sizes=="64" & 
#     results_cor_usemean_varmean_varsisd$i == "i == 2.67", ], 
#   aes(var_mean, var_si_sd)) + 
#   geom_point() + 
#   xlim(0, 26) +
#   ylim(0, 26) +
#   geom_abline(slope=1, intercept=0) + 
#   # facet_wrap(~i, labeller = label_parsed) + 
#   xlab("variance of mean") + 
#   ylab("variance of i*sd") +
#   theme_minimal_grid(font_size=8) + 
#   theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), 
#         plot.title=element_text(size=8)) + 
#   ggtitle("number of causal loci: 64")
# # dev.off()
# # pdf("view_usefulness/plots/test1.pdf")
# p7 <-
# ggplot(results_cor_usemean_varmean_varsisd[
#   results_cor_usemean_varmean_varsisd$effective_marker_sizes=="1024" & 
#     results_cor_usemean_varmean_varsisd$i == "i == 2.67", ], 
#   aes(var_mean, var_si_sd)) + 
#   geom_point() + 
#   xlim(0, 260) +
#   ylim(0, 260) +
#   geom_abline(slope=1, intercept=0) + 
#   # facet_wrap(~i, labeller = label_parsed) + 
#   xlab("variance of mean") + 
#   ylab("variance of i*sd") +
#   theme_minimal_grid(font_size=8) + 
#   theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), 
#         plot.title=element_text(size=8)) + 
#   ggtitle("number of causal loci: 1024")
# # dev.off()

# prow <- plot_grid(p5, p6, p7, nrow=1, labels=c("b", "c", "d"))

save_plot("view_usefulness/plots/cor_fam_use_mean_var_sisd_mean.pdf", 
          # prow,
          plot_grid(p4
                    # , p5, ncol=1, labels="auto", label_size = 7, align="vh", axis="lr"
                    ),
          base_width=3.5, base_height=2.15)
# save_plot("view_usefulness/plots/t_value.pdf", 
#           # prow,
#           plot_grid(p5),
#           base_width=3.5, base_height=2.15)


# predict fam mean bv
# predict fam sd bv with two methods 
# predict fam use with/without sd 
# corre between mean and sd*si true bv with slope 1 



# following is archived
{
# predict fam mean bv
cor_fammean_RR <- vector(mode="list", length=length(h2s))
cor_fammean_RR <- lapply(cor_fammean_RR, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (h in 1:20){
      cor_fammean_RR[[i]][[j]][h] = 
        cor(BV_fammean[[i]][[j]][, h], predY_RR_fammean[[i]][[j]][, h])
    }
  }
}
cor_fammean_RR_df <- data.frame(h2s = rep(h2s, each=(6*20)), 
                                effective_marker_sizes = 
                                  as.character(rep(rep(effective_marker_sizes, each=20), 5)), 
                                cor = rep(NA, (6*5*20)))
# cor_fammean_RR_df$h2s <- factor(cor_fammean_RR_df$h2s, levels=as.factor(h2s))
cor_fammean_RR_df$effective_marker_sizes <- 
  factor(cor_fammean_RR_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    cor_fammean_RR_df$cor[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] = 
      cor_fammean_RR[[i]][[j]]
  }
}
pdf(paste("view_usefulness/plots/", "fammean_RR_cor.pdf", sep=""))
ggplot(cor_fammean_RR_df, aes(as.numeric(effective_marker_sizes), cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~h2s) + 
  xlab("number of effective QTL") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  ggtitle("accuracy of predicting family mean of BV, rrBLUP")
dev.off()
pdf(paste("view_usefulness/plots/", "fammean_RR_cor_xaxis_h2s.pdf", sep=""))
ggplot(cor_fammean_RR_df, aes(h2s, cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~effective_marker_sizes) + 
  xlab("heritability") + 
  ylab("accuracy") + 
  ggtitle("accuracy of predicting family mean of BV, rrBLUP")
dev.off()

pdf(paste("view_usefulness/plots/", "fammean_RR_cor_2.pdf", sep=""), width=10)
ggplot(cor_fammean_RR_df[cor_fammean_RR_df$h2s %in% c(0.9, 0.5, 0.1), ], 
       aes(as.numeric(effective_marker_sizes), cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~h2s) + 
  xlab("number of effective QTL") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  ggtitle("accuracy of predicting family mean of BV, rrBLUP")
dev.off()




cor_fammean <- vector(mode="list", length=length(h2s))
cor_fammean <- lapply(cor_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (h in 1:20){
      cor_fammean[[i]][[j]][h] = 
        cor(BV_fammean[[i]][[j]][, h], predY_fammean[[i]][[j]][, h])
    }
  }
}
cor_fammean_df <- data.frame(h2s = rep(paste("h^2 == ", h2s, sep=""), each=(6*20)), 
                                effective_marker_sizes = 
                                  as.character(rep(rep(effective_marker_sizes, each=20), 5)), 
                                cor = rep(NA, (6*5*20)))
cor_fammean_df$effective_marker_sizes <-
  factor(cor_fammean_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# cor_fammean_df$h2s <- factor(cor_fammean_df$h2s, levels = paste("h2 = ", h2s, sep=""))
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    cor_fammean_df$cor[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] =
      cor_fammean[[i]][[j]]
  }
}

# saveRDS(cor_fammean_df, "view_usefulness/cor_fammean_df.rds")
# saveRDS(lbs, "view_usefulness/lbs.rds")
pdf(paste("view_usefulness/plots/", "fammean_cor.pdf", sep=""))
ggplot(cor_fammean_df, aes(as.numeric(effective_marker_sizes), cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~h2s, labeller = as_labeller(lbs, label_parsed)) + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  scale_fill_discrete(labels=parse(text=lbs)) + 
  # theme(strip.text = ggtext::element_markdown())
  ggtitle("accuracy of predicting family mean of BV, BayesC")
dev.off()
pdf(paste("view_usefulness/plots/", "fammean_cor_xaxis_h2s.pdf", sep=""))
ggplot(cor_fammean_df, aes(h2s, cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~effective_marker_sizes) + 
  xlab("heritability") + 
  ylab("accuracy") + 
  ggtitle("accuracy of predicting family mean of BV, BayesC")
dev.off()



# predict fam sd bv with predicted fam sd bv
# predict fam sd bv with number of heterozygous markers of parents 
cor_famsd_heter_markers <- vector(mode="list", length=length(h2s))
cor_famsd_heter_markers <- lapply(cor_famsd_heter_markers, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (h in 1:20){
      cor_famsd_heter_markers[[i]][[j]][h] = 
        cor(sqrt(BV_famvar[[i]][[j]][, h]), parents_heter_sites_noncausal[[j]][, h])
    }
  }
}

cor_famsd_RR <- vector(mode="list", length=length(h2s))
cor_famsd_RR <- lapply(cor_famsd_RR, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (h in 1:20){
      cor_famsd_RR[[i]][[j]][h] = 
        cor(sqrt(BV_famvar[[i]][[j]][, h]), sqrt(predY_RR_famvar[[i]][[j]][, h]))
    }
  }
}

cor_famsd_RR_df <- data.frame(h2s = rep(h2s, each=(6*20)),
                              effective_marker_sizes =
                                as.character(rep(rep(effective_marker_sizes, each=20), 5)),
                              cor = rep(NA, (6*5*20)),
                              cor_het_mar = rep(NA, (6*5*20)))
cor_famsd_RR_df$effective_marker_sizes <-
  factor(cor_famsd_RR_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    cor_famsd_RR_df$cor[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] =
      cor_famsd_RR[[i]][[j]]
    cor_famsd_RR_df$cor_het_mar[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] =
      cor_famsd_heter_markers[[i]][[j]]
  }
}
pdf(paste("view_usefulness/plots/", "famsd_RR_cor.pdf", sep=""))
ggplot(cor_famsd_RR_df, aes(as.numeric(effective_marker_sizes), y=cor)) +
  geom_point() +
  geom_smooth(method="loess") +
  facet_wrap(~h2s) +
  xlab("number of causal loci") +
  ylab("accuracy") +
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  ggtitle("accuracy of predicting family SD of BV, rrBLUP")
dev.off()
pdf(paste("view_usefulness/plots/", "famsd_RR_cor_xaxis_h2s.pdf", sep=""))
ggplot(cor_famsd_RR_df, aes(h2s, y=cor)) +
  geom_point() +
  geom_smooth(method="loess") +
  facet_wrap(~effective_marker_sizes) +
  xlab("heritability") +
  ylab("accuracy") +
  ggtitle("accuracy of predicting family SD of BV, rrBLUP")
dev.off()
pdf(paste("view_usefulness/plots/", "famsd_RR_cor_two_way.pdf", sep=""))
ggplot(cor_famsd_RR_df, aes(as.numeric(effective_marker_sizes))) +
  geom_point(aes(y=cor, colour="prediction from model")) +
  geom_point(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents")) +
  geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
  geom_smooth(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents"), 
              method="loess") +
  facet_wrap(~h2s) +
  theme(legend.position="bottom") +
  xlab("number of causal loci") +
  ylab("accuracy/correlation") +
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  ggtitle("accuracy of predicting family SD of BV, rrBLUP")
dev.off()
pdf(paste("view_usefulness/plots/", "famsd_RR_cor_xaxis_h2s_two_way.pdf", sep=""))
ggplot(cor_famsd_RR_df, aes(h2s)) +
  geom_point(aes(y=cor, colour="prediction from model")) +
  geom_point(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents")) +
  geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
  geom_smooth(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents"), 
              method="loess") +
  facet_wrap(~effective_marker_sizes) +
  theme(legend.position="bottom") +
  xlab("heritability") +
  ylab("accuracy/correlation") +
  ggtitle("accuracy of predicting family SD of BV, rrBLUP")
dev.off()

pdf(paste("view_usefulness/plots/", "famsd_RR_cor_2.pdf", sep=""), width=10)
ggplot(cor_famsd_RR_df[cor_famsd_RR_df$h2s %in% c(0.9, 0.5, 0.1), ], 
       aes(as.numeric(effective_marker_sizes), y=cor)) +
  geom_point() +
  geom_smooth(method="loess") +
  facet_wrap(~h2s) +
  xlab("number of causal loci") +
  ylab("accuracy") +
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  ggtitle("accuracy of predicting family SD of BV, rrBLUP")
dev.off()
pdf(paste("view_usefulness/plots/", "famsd_RR_cor_two_way_2.pdf", sep=""), width=10)
ggplot(cor_famsd_RR_df[cor_famsd_RR_df$h2s %in% c(0.9, 0.5, 0.1), ], 
       aes(as.numeric(effective_marker_sizes))) +
  geom_point(aes(y=cor, colour="prediction from model")) +
  geom_point(aes(y=cor_het_mar, colour="correlation with number of hetero markers of parents")) +
  geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
  geom_smooth(aes(y=cor_het_mar, colour="correlation with number of hetero markers of parents"), 
              method="loess") +
  facet_wrap(~h2s) +
  theme(legend.position="bottom") +
  xlab("number of causal loci") +
  ylab("accuracy/correlation") +
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  ggtitle("accuracy of predicting family SD of BV, rrBLUP")
dev.off()



cor_famsd <- vector(mode="list", length=length(h2s))
cor_famsd <- lapply(cor_famsd, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (h in 1:20){
      cor_famsd[[i]][[j]][h] = 
        cor(sqrt(BV_famvar[[i]][[j]][, h]), sqrt(predY_famvar[[i]][[j]][, h]))
    }
  }
}

cor_famsd_df <- data.frame(h2s = rep(paste("h^2 == ", h2s, sep=""), each=(6*20)),
                           effective_marker_sizes =
                             as.character(rep(rep(effective_marker_sizes, each=20), 5)),
                           cor = rep(NA, (6*5*20)),
                           cor_het_mar = rep(NA, (6*5*20)))
cor_famsd_df$effective_marker_sizes <-
  factor(cor_famsd_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# cor_famsd_df$effective_marker_sizes <- as.numeric(cor_famsd_df$effective_marker_sizes)
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    cor_famsd_df$cor[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] =
      cor_famsd[[i]][[j]]
    cor_famsd_df$cor_het_mar[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] =
      cor_famsd_heter_markers[[i]][[j]]
  }
}
pdf(paste("view_usefulness/plots/", "famsd_cor.pdf", sep=""))
ggplot(cor_famsd_df, aes(as.numeric(effective_marker_sizes), y=cor)) +
  geom_point() +
  geom_smooth(method="loess") +
  facet_wrap(~h2s) +
  xlab("number of causal loci") +
  ylab("accuracy") +
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  ggtitle("accuracy of predicting family SD of BV, BayesC")
dev.off()
pdf(paste("view_usefulness/plots/", "famsd_cor_xaxis_h2s.pdf", sep=""))
ggplot(cor_famsd_df, aes(h2s, y=cor)) +
  geom_point() +
  geom_smooth(method="loess") +
  facet_wrap(~effective_marker_sizes) +
  xlab("heritability") +
  ylab("accuracy") +
  ggtitle("accuracy of predicting family SD of BV, BayesC")
dev.off()
pdf(paste("view_usefulness/plots/", "famsd_cor_two_way.pdf", sep=""))
ggplot(cor_famsd_df, aes(as.numeric(effective_marker_sizes))) +
  geom_point(aes(y=cor, colour="prediction from model")) +
  geom_point(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents")) +
  geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
  geom_smooth(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents"), 
              method="loess") +
  facet_wrap(~h2s) +
  theme(legend.position="bottom") +
  xlab("number of causal loci") +
  ylab("accuracy/correlation") +
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  ggtitle("accuracy of predicting family SD of BV, BayesC")
dev.off()
pdf(paste("view_usefulness/plots/", "famsd_cor_xaxis_h2s_two_way.pdf", sep=""))
ggplot(cor_famsd_df, aes(h2s)) +
  geom_point(aes(y=cor, colour="prediction from model")) +
  geom_point(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents")) +
  geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
  geom_smooth(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents"), 
              method="loess") +
  facet_wrap(~effective_marker_sizes) +
  theme(legend.position="bottom") +
  xlab("heritability") +
  ylab("accuracy/correlation") +
  ggtitle("accuracy of predicting family SD of BV, BayesC")
dev.off()



# cor_fammean_df$cor_het_mar <- cor_famsd_df$cor_het_mar
# cor_fammean_df$type <- "mean"
# cor_famsd_df$type <- "sd"


# lbs <- setNames(
#   paste("''*h^2*' = ", h2s, "'", sep=""), levels(cor_fammean_df$h2s)
# )
p1 <- ggplot(cor_fammean_df[cor_fammean_df$h2s %in% c("h^2 == 0.9", "h^2 == 0.5", "h^2 == 0.1"), ], 
             aes(as.numeric(effective_marker_sizes), cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~h2s, labeller = label_parsed) + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  # ggtitle("accuracy of predicting family mean of BV, BayesC") + 
  theme_minimal_grid(font_size=10) +
  theme(axis.title.x=element_blank(),
       axis.text.x=element_blank(),
       axis.ticks.x=element_blank())

# save_plot(paste("view_usefulness/plots/", "fammean_cor_2.pdf", sep=""),
#           p1, 
#           base_width=6.5, base_height=2.17)

p2 <- ggplot(cor_famsd_df[cor_famsd_df$h2s %in% c("h^2 == 0.9", "h^2 == 0.5", "h^2 == 0.1"), ], 
               aes(as.numeric(effective_marker_sizes))) +
  geom_point(aes(y=cor, colour="prediction from model")) +
  geom_point(aes(y=cor_het_mar, colour="prediction from parental heterozygosity")) +
  geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
  geom_smooth(aes(y=cor_het_mar, colour="prediction from parental heterozygosity"), 
              method="loess") +
  facet_wrap(~h2s, labeller = label_parsed) + 
  # ylim(-0.25, 1) +
  xlab("number of causal loci") +
  ylab("accuracy") +
  ylim(0, 1) + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  theme_minimal_grid(font_size=10) +
  theme(strip.text.x = element_blank(), 
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=1))
  #       axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())
  # theme(strip.background = element_blank())

prow <- plot_grid(p1, 
                  p2, 
                  labels="auto", ncol=1, rel_heights = c(0.85, 1))

save_plot(paste("view_usefulness/plots/", "fammean_famsd_cor.pdf", sep=""), 
          plot_grid(prow), 
          base_width=6.5, base_height = 4.33)
}

# following is archived
{
# predict fam use with/without sd 
# cor(BV_famuse, predY_RR_fammean)
cor_famuse_RR_fammean <- vector(mode="list", length=length(si))
cor_famuse_RR_fammean <- lapply(cor_famuse_RR_fammean, FUN=function(x){
  vector(mode="list", length=length(h2s))
})
for (h in 1:length(si)){
  cor_famuse_RR_fammean[[h]] = lapply(cor_famuse_RR_fammean[[h]], FUN=function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
}
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        cor_famuse_RR_fammean[[h]][[i]][[j]][k] = 
          cor(BV_famuse[[h]][[i]][[j]][, k], predY_RR_fammean[[i]][[j]][, k])
      }
    }
  }
}
# cor(BV_famuse, predY_RR_famuse)
cor_famuse_RR_famuse <- vector(mode="list", length=length(si))
cor_famuse_RR_famuse <- lapply(cor_famuse_RR_famuse, FUN=function(x){
  vector(mode="list", length=length(h2s))
})
for (h in 1:length(si)){
  cor_famuse_RR_famuse[[h]] = lapply(cor_famuse_RR_famuse[[h]], FUN=function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
}
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        cor_famuse_RR_famuse[[h]][[i]][[j]][k] = 
          cor(BV_famuse[[h]][[i]][[j]][, k], predY_RR_famuse[[h]][[i]][[j]][, k])
      }
    }
  }
}
cor_famuse_RR_df <- data.frame(si = rep(si, each=(5*6*20)), 
                               h2s = rep(rep(h2s, each=(6*20)), 4), 
                               effective_marker_sizes = 
                                 as.character(rep(rep(effective_marker_sizes, each=20), (4*5))), 
                               cor_mean = rep(NA, (4*6*5*20)),
                               cor_use = rep(NA, (4*6*5*20)))
cor_famuse_RR_df$effective_marker_sizes <- 
  factor(cor_famuse_RR_df$effective_marker_sizes, levels=effective_marker_sizes)
cor_famuse_RR_df$si <- round(cor_famuse_RR_df$si, 2)
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      cor_famuse_RR_df$cor_mean[((h-1)*5*6*20 + (i-1)*6*20 + (j-1)*20 + 1) : 
                                  ((h-1)*5*6*20 + (i-1)*6*20 + (j-1)*20 + 20)] = 
        cor_famuse_RR_fammean[[h]][[i]][[j]]
      cor_famuse_RR_df$cor_use[((h-1)*5*6*20 + (i-1)*6*20 + (j-1)*20 + 1) : 
                                  ((h-1)*5*6*20 + (i-1)*6*20 + (j-1)*20 + 20)] = 
        cor_famuse_RR_famuse[[h]][[i]][[j]]
    }
  }
}
pdf(paste("view_usefulness/plots/", "famuse_RR_cor.pdf", sep=""), width=14, height=10)
ggplot(cor_famuse_RR_df, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=cor_mean, colour="predicted from family mean BV")) + 
  geom_point(aes(y=cor_use, colour="predicted from family usefulness BV")) + 
  geom_smooth(aes(y=cor_mean, colour="predicted from family mean BV"), method="loess") + 
  geom_smooth(aes(y=cor_use, colour="predicted from family usefulness BV"), method="loess") + 
  facet_grid(si ~ h2s) + 
  theme(legend.position="bottom") + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  ggtitle("accuracy of predicting family usefulness of BV, rrBLUP")
dev.off()
pdf(paste("view_usefulness/plots/", "famuse_RR_cor_xaxis_h2s.pdf", sep=""), 
    width=14, height=10)
ggplot(cor_famuse_RR_df, aes(h2s)) + 
  geom_point(aes(y=cor_mean, colour="predicted from family mean BV")) + 
  geom_point(aes(y=cor_use, colour="predicted from family usefulness BV")) + 
  geom_smooth(aes(y=cor_mean, colour="predicted from family mean BV"), method="loess") + 
  geom_smooth(aes(y=cor_use, colour="predicted from family usefulness BV"), method="loess") + 
  facet_grid(si~effective_marker_sizes) + 
  theme(legend.position="bottom") + 
  xlab("heritability") + 
  ylab("accuracy") + 
  ggtitle("accuracy of predicting family usefulness of BV, rrBLUP")
dev.off()
pdf(paste("view_usefulness/plots/", "famuse_RR_cor_2.pdf", sep=""), width=14, height=10)
ggplot(cor_famuse_RR_df[cor_famuse_RR_df$h2s %in% c(0.9, 0.5, 0.1) & 
                          cor_famuse_RR_df$si %in% c(1.4, 2.42), ], 
       aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=cor_mean, colour="predicted from family mean BV")) + 
  geom_point(aes(y=cor_use, colour="predicted from family usefulness BV")) + 
  geom_smooth(aes(y=cor_mean, colour="predicted from family mean BV"), method="loess") + 
  geom_smooth(aes(y=cor_use, colour="predicted from family usefulness BV"), method="loess") + 
  facet_grid(si ~ h2s) + 
  theme(legend.position="bottom") + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  ggtitle("accuracy of predicting family usefulness of BV, rrBLUP")
dev.off()

# cor(BV_famuse, predY_fammean)
cor_famuse_fammean <- vector(mode="list", length=length(si))
cor_famuse_fammean <- lapply(cor_famuse_fammean, FUN=function(x){
  vector(mode="list", length=length(h2s))
})
for (h in 1:length(si)){
  cor_famuse_fammean[[h]] = lapply(cor_famuse_fammean[[h]], FUN=function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
}
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        cor_famuse_fammean[[h]][[i]][[j]][k] = 
          cor(BV_famuse[[h]][[i]][[j]][, k], predY_fammean[[i]][[j]][, k])
      }
    }
  }
}
# cor(BV_famuse, predY_famuse)
cor_famuse_famuse <- vector(mode="list", length=length(si))
cor_famuse_famuse <- lapply(cor_famuse_famuse, FUN=function(x){
  vector(mode="list", length=length(h2s))
})
for (h in 1:length(si)){
  cor_famuse_famuse[[h]] = lapply(cor_famuse_famuse[[h]], FUN=function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
}
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        cor_famuse_famuse[[h]][[i]][[j]][k] = 
          cor(BV_famuse[[h]][[i]][[j]][, k], predY_famuse[[h]][[i]][[j]][, k])
      }
    }
  }
}
cor_famuse_df <- data.frame(si = rep(paste("i == ", round(si, 2), sep=""), each=(5*6*20)), 
                               h2s = rep(rep(paste("h^2 == ", h2s, sep=""), each=(6*20)), 4), 
                               effective_marker_sizes = 
                                 as.character(rep(rep(effective_marker_sizes, each=20), (4*5))), 
                               cor_mean = rep(NA, (4*6*5*20)),
                               cor_use = rep(NA, (4*6*5*20)))
cor_famuse_df$effective_marker_sizes <- 
  factor(cor_famuse_df$effective_marker_sizes, levels=effective_marker_sizes)
# cor_famuse_df$si <- round(cor_famuse_df$si, 2)

for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      cor_famuse_df$cor_mean[((h-1)*5*6*20 + (i-1)*6*20 + (j-1)*20 + 1) : 
                                  ((h-1)*5*6*20 + (i-1)*6*20 + (j-1)*20 + 20)] = 
        cor_famuse_fammean[[h]][[i]][[j]]
      cor_famuse_df$cor_use[((h-1)*5*6*20 + (i-1)*6*20 + (j-1)*20 + 1) : 
                                 ((h-1)*5*6*20 + (i-1)*6*20 + (j-1)*20 + 20)] = 
        cor_famuse_famuse[[h]][[i]][[j]]
    }
  }
}
{pdf(paste("view_usefulness/plots/", "famuse_cor.pdf", sep=""), width=14, height=10)
ggplot(cor_famuse_df, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=cor_mean, colour="predicted from family mean BV")) + 
  geom_point(aes(y=cor_use, colour="predicted from family usefulness BV")) + 
  geom_smooth(aes(y=cor_mean, colour="predicted from family mean BV"), method="loess") + 
  geom_smooth(aes(y=cor_use, colour="predicted from family usefulness BV"), method="loess") + 
  facet_grid(si~h2s) + 
  theme(legend.position="bottom") + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  ggtitle("accuracy of predicting family usefulness of BV, BayesC")
dev.off()
pdf(paste("view_usefulness/plots/", "famuse_cor_xaxis_h2s.pdf", sep=""), width=14, height=10)
ggplot(cor_famuse_df, aes(h2s)) + 
  geom_point(aes(y=cor_mean, colour="predicted from family mean BV")) + 
  geom_point(aes(y=cor_use, colour="predicted from family usefulness BV")) + 
  geom_smooth(aes(y=cor_mean, colour="predicted from family mean BV"), method="loess") + 
  geom_smooth(aes(y=cor_use, colour="predicted from family usefulness BV"), method="loess") + 
  facet_grid(si~effective_marker_sizes) + 
  theme(legend.position="bottom") + 
  xlab("heritability") + 
  ylab("accuracy") + 
  ggtitle("accuracy of predicting family usefulness of BV, BayesC")
dev.off()
pdf(paste("view_usefulness/plots/", "famuse_cor_2.pdf", sep=""), width=14, height=10)
ggplot(cor_famuse_df[cor_famuse_df$h2s %in% c("h2 = 0.9", "h2 = 0.5", "h2 = 0.1") & 
                       cor_famuse_df$si %in% c("si = 1.4", "si = 2.42"), ], 
       aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=cor_mean, colour="predicted from family mean BV")) + 
  geom_point(aes(y=cor_use, colour="predicted from family usefulness BV")) + 
  geom_smooth(aes(y=cor_mean, colour="predicted from family mean BV"), method="loess") + 
  geom_smooth(aes(y=cor_use, colour="predicted from family usefulness BV"), method="loess") + 
  facet_grid(si~h2s) + 
  theme(legend.position="bottom") + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  ggtitle("accuracy of predicting family usefulness of BV, BayesC")
dev.off()
}

p3 <- ggplot(cor_famuse_df[cor_famuse_df$h2s %in% c("h^2 == 0.9", "h^2 == 0.5", "h^2 == 0.1") & 
                             cor_famuse_df$si %in% c("i == 1.4", "i == 2.42"), ], 
             aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=cor_mean, colour="predicted from family mean BV")) + 
  geom_point(aes(y=cor_use, colour="predicted from family usefulness BV")) + 
  geom_smooth(aes(y=cor_mean, colour="predicted from family mean BV"), method="loess") + 
  geom_smooth(aes(y=cor_use, colour="predicted from family usefulness BV"), method="loess") + 
  facet_grid(si~h2s, labeller = label_parsed) +
  # facet_wrap(~h2s, labeller = as_labeller(lbs, label_parsed)) + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  guides(color=guide_legend(nrow=1)) + 
  # ggtitle("accuracy of predicting family usefulness of BV, BayesC") + 
  theme_minimal_grid(font_size=10) +
  theme(legend.position="bottom") 
save_plot(paste("view_usefulness/plots/", "famuse_cor_2.pdf", sep=""), 
          plot_grid(p3), 
          base_width=6.5, base_height=4.33)
}



# following is archived
{# corre between usefulness and mean
# cor(BV_famuse, BV_fammean)
cor_BV_famuse_fammean <- vector(mode="list", length=length(si))
cor_BV_famuse_fammean <- lapply(cor_BV_famuse_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(si)){
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      cor_BV_famuse_fammean[[i]][[j]][k] = 
        cor(BV_famuse[[i]][[1]][[j]][, k], BV_fammean[[1]][[j]][, k])
    }
  }
}
cor_BV_famuse_fammean_df <- data.frame(si=rep(si, each=(6*20)), 
                                       effective_marker_sizes=
                                         rep(rep(effective_marker_sizes, each=20), 4), 
                                       cor=rep(NA, 4*6*20))
cor_BV_famuse_fammean_df$effective_marker_sizes <- 
  factor(cor_BV_famuse_fammean_df$effective_marker_sizes, levels=effective_marker_sizes)
cor_BV_famuse_fammean_df$si <- round(cor_BV_famuse_fammean_df$si, 2)
for (i in 1:length(si)){
  for (j in 1:length(effective_marker_sizes)){
    cor_BV_famuse_fammean_df$cor[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] = 
      cor_BV_famuse_fammean[[i]][[j]]
  }
}
pdf(paste("view_usefulness/plots/", "BV_famuse_fammean.pdf", sep=""))
ggplot(cor_BV_famuse_fammean_df, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=cor)) + 
  geom_smooth(aes(y=cor), method="loess") + 
  facet_wrap(~si) + 
  theme(legend.position="bottom") + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  theme_minimal_grid(font_size=10)
  # ggtitle("correlation between BV family usefulness and BV family mean")
dev.off()




fam_mean_sd <- data.frame(h2s = rep(h2s, each=(6*length(file_names))), 
                          effective_marker_sizes = 
                            rep(rep(effective_marker_sizes, each=length(file_names)), 5), 
                          family = rep(file_names, (6*5)), 
                          matrix(NA, nrow=(30*length(file_names)), ncol=120))
fam_mean_sd$effective_marker_sizes <- 
  factor(fam_mean_sd$effective_marker_sizes, levels=as.factor(effective_marker_sizes))

phenotype_types <- c()
for (i in 1:20){
  phenotype_types[i] = paste("BV_fammean_", i, sep="")
}
for (i in 21:40){
  phenotype_types[i] = paste("BV_famsd_", (i-20), sep="")
}
for (i in 41:60){
  phenotype_types[i] = paste("predY_RR_fammean_", (i-40), sep="")
}
for (i in 61:80){
  phenotype_types[i] = paste("predY_RR_famsd_", (i-60), sep="")
}
for (i in 81:100){
  phenotype_types[i] = paste("predY_fammean_", (i-80), sep="")
}
for (i in 101:120){
  phenotype_types[i] = paste("predY_famsd_", (i-100), sep="")
}
names(fam_mean_sd)[4:123] <- phenotype_types

for (i in 1:length(BV_fammean)){
  for (j in 1:length(BV_fammean[[i]])){
    fam_mean_sd[((i-1)*6*length(file_names) + (j-1)*length(file_names) + 1) : 
                  ((i-1)*6*length(file_names) + (j-1)*length(file_names) + 520), 4:43] = 
      cbind(BV_fammean[[i]][[j]], sqrt(BV_famvar[[i]][[j]]))
  }
}
for (i in 1:length(predY_RR_fammean)){
  for (j in 1:length(predY_RR_fammean[[i]])){
    fam_mean_sd[((i-1)*6*length(file_names) + (j-1)*length(file_names) + 1) : 
                  ((i-1)*6*length(file_names) + (j-1)*length(file_names) + 520), 44:83] = 
      cbind(predY_RR_fammean[[i]][[j]], sqrt(predY_RR_famvar[[i]][[j]]))
  }
}
for (i in 1:length(predY_fammean)){
  for (j in 1:length(predY_fammean[[i]])){
    fam_mean_sd[((i-1)*6*length(file_names) + (j-1)*length(file_names) + 1) : 
                  ((i-1)*6*length(file_names) + (j-1)*length(file_names) + 520), 84:123] = 
      cbind(predY_fammean[[i]][[j]], sqrt(predY_famvar[[i]][[j]]))
  }
}

saveRDS(fam_mean_sd, "view_usefulness/fam_mean_sd.rds")
fam_mean_sd <- readRDS("view_usefulness/fam_mean_sd.rds")


fam_mean_sd <- fam_mean_sd[, 1:43]
fam_mean_sd <- fam_mean_sd[fam_mean_sd$h2s==0.9, ][, 2:43]
fam_mean_sd$effective_marker_sizes <- as.numeric(as.character(fam_mean_sd$effective_marker_sizes))
fam_mean_sd$family <- as.character(fam_mean_sd$family)

fam_mean_sd_all <- data.frame(h2s=NA, number_of_QTL=NA, family=NA, replic=NA, mean=NA, sd=NA)
for (i in 1:nrow(fam_mean_sd)){
  fam_mean_sd_all = rbind(fam_mean_sd_all, 
                          data.frame(h2s=rep(NA, 20),
                                     number_of_QTL=rep(fam_mean_sd$effective_marker_sizes[i], 20), 
                                     family=rep(fam_mean_sd$family[i], 20), 
                                     replic=1:20, 
                                     mean=as.numeric(fam_mean_sd[i, 3:22][1, ]), 
                                     sd=as.numeric(fam_mean_sd[i, 23:42][1, ])))
}
fam_mean_sd_all <- fam_mean_sd_all[-1, -1]
fam_mean_sd_all[, c("parent1", "parent2")] <- 
  t(sapply(fam_mean_sd_all$family[530:535], FUN=function(x){
    strsplit(x, "_")[[1]]
  }))
fam_mean_sd_all$parents_type <- "random"
fam_mean_sd_all <- fam_mean_sd_all[, c(1, 3, 6, 7, 8, 4, 5)]

fam_mean_sd_use <- fam_mean_sd_all[rep(1:nrow(fam_mean_sd_all), 4), ]
fam_mean_sd_use$si <- rep(si, each=nrow(fam_mean_sd_all))
fam_mean_sd_use$si_sd <- fam_mean_sd_use$si * fam_mean_sd_use$sd
fam_mean_sd_use$use <- fam_mean_sd_use$mean + fam_mean_sd_use$si * fam_mean_sd_use$sd

fam_cor <- data.frame(number_of_QTL=rep(as.character(effective_marker_sizes), each=4*20), 
                      si=rep(rep(paste("i = ", round(si, 2), sep=""), 6), each=20),
                      replic=rep(1:20, 4*6),
                      cor_use_mean=rep(NA, 4*6*20))

for (i in 1:length(effective_marker_sizes)){
  for (j in 1:length(si)){
    for (k in 1:20){
      temp = fam_mean_sd_use[(fam_mean_sd_use$number_of_QTL == effective_marker_sizes[i]) & 
                               (fam_mean_sd_use$si == si[j]) & 
                               (fam_mean_sd_use$replic == k), ]
      fam_cor[(i-1)*4*20 + (j-1)*20 + k, "cor_use_mean"] = cor(temp$mean, temp$use)
    }
  }
}
fam_cor$number_of_QTL <- factor(fam_cor$number_of_QTL, levels=as.factor(effective_marker_sizes))
# fam_cor$si <- factor(round(fam_cor$si, 2), levels=round(si, 2))
fam_cor$replic <- factor(fam_cor$replic, levels=1:20)



fam_var_mean_var_sd <- aggregate(fam_mean_sd_use$mean, list(fam_mean_sd_use$number_of_QTL, 
                                                            fam_mean_sd_use$replic, 
                                                            fam_mean_sd_use$si), var)
colnames(fam_var_mean_var_sd) <- c("number_of_QTL", "replic", "si", "var_mean")
fam_var_mean_var_sd$var_si_sd <- aggregate(fam_mean_sd_use$si_sd, list(fam_mean_sd_use$number_of_QTL, 
                                                                       fam_mean_sd_use$replic, 
                                                                       fam_mean_sd_use$si), var)[, 4]
fam_var_mean_var_sd$si <- round(fam_var_mean_var_sd$si, 2)
fam_var_mean_var_sd$si <- paste("i = ", (fam_var_mean_var_sd$si), sep="")

fam_var_mean_var_sd_1024 <- fam_var_mean_var_sd[fam_var_mean_var_sd$number_of_QTL==1024, ]
summary(fam_var_mean_var_sd_1024$var_mean)
seqnc <- seq(100, 255, length.out=5)
fam_var_mean_var_sd_1024$group <- .bincode(fam_var_mean_var_sd_1024$var_mean, seqnc)
xlab <- cut(seq(100, 255, length.out=5), seqnc)
xlab <- as.character(xlab)
fam_var_mean_var_sd_1024$xlab <- NA
for (i in 1:nrow(fam_var_mean_var_sd_1024)){
  fam_var_mean_var_sd_1024$xlab[i] = xlab[fam_var_mean_var_sd_1024$group[i]+1]
}
fam_var_mean_var_sd_1024$xlab <- factor(fam_var_mean_var_sd_1024$xlab, levels=xlab)

p4 <- ggplot(fam_cor, aes(x=as.numeric(number_of_QTL), y=cor_use_mean)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~si) + 
  xlab("number of causal loci") + 
  ylab("correlation between BV mean and usefulness") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  theme_minimal_grid(font_size=8)
p5 <- ggplot(fam_var_mean_var_sd_1024, aes(var_mean, var_si_sd)) + 
  geom_point() + 
  facet_wrap(~si) + 
  xlab("variance of mean") + 
  ylab("variance of i*sd") +
  theme_minimal_grid(font_size=8)

save_plot("view_usefulness/plots/cor_fam_use_mean_var_sisd_mean.pdf", 
          plot_grid(p4, p5, nrow=1, labels="auto"),
          base_width=6.5, base_height = 4.33)
}

