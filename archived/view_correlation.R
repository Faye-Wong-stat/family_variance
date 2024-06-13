setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)
sf <- c(0.8, 0.99)
b <- qnorm(sf, 0, 1)
si <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)
select_number <- c(200, 100, 50, 25, 5)

file_names_results <- list.files("simulate_phenotypes_crosses/", pattern="_result_df.rds")
file_names <- gsub("_result_df.rds", "", file_names_results)
parents_names <- matrix(NA, nrow=length(file_names), ncol=2)
parents_names[, 1] <- gsub("_.+", "", file_names)
parents_names[, 2] <- gsub(".+_", "", file_names)
rownames(parents_names) <- file_names

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



results <- readRDS("view_usefulness/results.rds")
file_names_results <- list.files("simulate_phenotypes_crosses2/", pattern="_result_df.rds")

temp <- lapply(file_names_results, FUN=function(x) {
  readRDS(paste("simulate_phenotypes_crosses2/", x, sep=""))
})
for (i in 1:length(temp)){
  temp[[i]]$BV_sd = sqrt(temp[[i]]$BV_var)
  temp[[i]]$predY_RR_sd = sqrt(temp[[i]]$predY_RR_var)
  temp[[i]]$predY_sd = sqrt(temp[[i]]$predY_var)
  temp[[i]] = temp[[i]][, -c(7, 9, 11)]
}

results$number_offspring <- 200
results <- results[, c(1, 11, 2:10)]

results_2 <- as.data.frame(matrix(NA, nrow=60*10000, ncol=11))
colnames(results_2) <- colnames(results)

for (i in 1:length(temp)){
  if(i %% 500 == 0){
    print(i)
  }
  results_2[((i-1)*60 + 1) : ((i-1)*60 + 60), ] = temp[[i]]
}
results_2$effective_marker_sizes <- factor(results_2$effective_marker_sizes, labels = effective_marker_sizes)

results <- rbind(results, results_2)
saveRDS(results, "view_correlation/results.rds")
results <- readRDS("view_correlation/results.rds")



# 900 = 5*5*3*20
results_cor <- matrix(NA, nrow=1500, ncol=9)
results_cor <- as.data.frame(results_cor)
colnames(results_cor) <- c("number_offspring", "h2s", "effective_marker_sizes", "trait_number", 
                           "mean_RR_cor", "mean_cor", "sd_RR_cor", "sd_cor", "sd_het_cor")
for (h in 1:length(select_number)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        # print((i-1)*5*20 + (j-1)*20 + k)
        A = results[results$number_offspring==select_number[h] & 
                      results$h2s == h2s[i] & 
                      results$effective_marker_sizes == effective_marker_sizes[j] & 
                      results$trait_number == k, ]
        results_cor[(h-1)*3*5*20 + (i-1)*5*20 + (j-1)*20 + k, ] = 
          c(select_number[h], h2s[i], effective_marker_sizes[j], k, 
            cor(A$BV_mean, A$predY_RR_mean), cor(A$BV_mean, A$predY_mean), 
            cor(A$BV_sd, A$predY_RR_sd), cor(A$BV_sd, A$predY_sd), 
            cor(A$BV_sd, parents_heter_sites_noncausal[[j]][, k]))
      }
    }
  }
}
results_cor$effective_marker_sizes <- factor(results_cor$effective_marker_sizes, 
                                             levels=as.factor(effective_marker_sizes))

# 5*5*3 = 75
results_cor_se <- matrix(NA, nrow=75, ncol=13)
results_cor_se <- as.data.frame(results_cor_se)
colnames(results_cor_se) <- c("number_offspring", "h2s", "effective_marker_sizes", 
                              "mean_cor_mean_RR", "mean_cor_se_RR", "sd_cor_mean_RR", "sd_cor_se_RR", 
                              "mean_cor_mean", "mean_cor_se", "sd_cor_mean", "sd_cor_se", 
                              "sd_het_cor_mean", "sd_het_cor_se")
for (h in 1:length(select_number)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      A = results_cor[results_cor$number_offspring==select_number[h] & 
                        results_cor$h2s == h2s[i] & 
                        results_cor$effective_marker_sizes == effective_marker_sizes[j], ]
      results_cor_se[(h-1)*3*5 + (i-1)*5 + j, ] = 
        c(select_number[h], h2s[i], effective_marker_sizes[j], 
          mean(A$mean_RR_cor), sd(A$mean_RR_cor)/sqrt(20), mean(A$sd_RR_cor), sd(A$sd_RR_cor)/sqrt(20), 
          mean(A$mean_cor), sd(A$mean_cor)/sqrt(20), mean(A$sd_cor), sd(A$sd_cor)/sqrt(20), 
          mean(A$sd_het_cor), sd(A$sd_het_cor)/sqrt(20))
    }
  }
}

p1 <- ggplot(results_cor_se, aes(number_offspring)) + 
  geom_point(aes(y=mean_cor_mean), color="blue") + 
  geom_errorbar(aes(ymin=mean_cor_mean-mean_cor_se, ymax=mean_cor_mean+mean_cor_se), 
                width=20, color="blue") + 
  geom_line(aes(y=mean_cor_mean), color="blue", linewidth=0.5) + 
  facet_grid(h2s~effective_marker_sizes) + 
  xlab("number of offspring used") + 
  ylab("prediction accuracy of family mean") + 
  ylim(0, 1) + 
  # ggtitle("accuracy of predicting family mean of BV, BayesC") + 
  theme_minimal_grid(font_size=10) #+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())
p2 <- ggplot(results_cor_se, aes(number_offspring)) + 
  geom_point(aes(y=sd_cor_mean, colour="prediction from model")) + 
  geom_errorbar(aes(ymin=sd_cor_mean-sd_cor_se, ymax=sd_cor_mean+sd_cor_se, 
                    colour="prediction from model"), width=20) + 
  geom_line(aes(y=sd_cor_mean, colour="prediction from model"), linewidth=0.5) + 
  geom_point(aes(y=sd_het_cor_mean, colour="prediction from parental heterozygosity")) + 
  geom_errorbar(aes(ymin=sd_het_cor_mean-sd_het_cor_se, ymax=sd_het_cor_mean+sd_het_cor_se, 
                    colour="prediction from parental heterozygosity"), width=20) + 
  geom_line(aes(y=sd_het_cor_mean, colour="prediction from parental heterozygosity"), linewidth=0.5) + 
  facet_grid(h2s~effective_marker_sizes) + 
  xlab("number of offspring used") + 
  ylab("prediction accuracy of family sd") + 
  ylim(0, 1) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=10) + 
  theme(legend.position = "bottom") +  
  guides(colour=guide_legend(title="prediction method", nrow=1))

save_plot(paste("view_correlation/plots/", "test1.pdf", sep=""), 
          # p1, 
          plot_grid(p1, p2, labels="auto", ncol=1, rel_heights = c(0.85, 1)),
          base_width=6.5, base_height = 8.66)



results_use <- results[rep(1:nrow(results), 2), ]
results_use$si <- rep(si, each=nrow(results))
results_use <- results_use[, c(1, 12, 2:11)]
results_use$BV_use <- results_use$BV_mean + results_use$si * results_use$BV_sd
results_use$predY_RR_use <- results_use$predY_RR_mean + results_use$si * results_use$predY_RR_sd
results_use$predY_use <- results_use$predY_mean + results_use$si * results_use$predY_sd

# 5*2*3*5*20 = 3000
results_cor_use <- matrix(NA, nrow=3000, ncol=10)
results_cor_use <- as.data.frame(results_cor_use)
colnames(results_cor_use) <- c("number_offspring", "si", "h2s", "effective_marker_sizes", "trait_number", 
                               "use_mean_RR_cor", "use_use_RR_cor", "use_mean_cor", "use_use_cor", 
                               "realuse_realmean_cor")
for (g in 1:length(select_number)){
  for (h in 1:length(si)){
    for (i in 1:length(h2s)){
      for (j in 1:length(effective_marker_sizes)){
        for (k in 1:20){
          A = results_use[results_use$number_offspring==select_number[g] & 
                            results_use$si == si[h] & 
                            results_use$h2s == h2s[i] & 
                            results_use$effective_marker_sizes == effective_marker_sizes[j] & 
                            results_use$trait_number == k, ]
          results_cor_use[(g-1)*2*3*5*20 + (h-1)*3*5*20 + (i-1)*5*20 + (j-1)*20 + k, ] = 
            c(select_number[g], si[h], h2s[i], effective_marker_sizes[j], k, 
              cor(A$BV_use, A$predY_RR_mean), cor(A$BV_use, A$predY_RR_use), 
              cor(A$BV_use, A$predY_mean), cor(A$BV_use, A$predY_use), 
              cor(A$BV_use, A$BV_mean))
        }
      }
    }
  }
}

results_cor_use$effective_marker_sizes <- factor(results_cor_use$effective_marker_sizes, 
                                                 levels=as.factor(effective_marker_sizes))
# results_cor_use_2 <- results_cor_use[results_cor_use$si==si[2], ]

# 5*2*3*5
results_cor_use_se <- matrix(NA, nrow=150, ncol=12)
results_cor_use_se <- as.data.frame(results_cor_use_se)
colnames(results_cor_use_se) <- c("number_offspring", "si", "h2s", "effective_marker_sizes", 
                                  "use_mean_RR_mean", "use_mean_RR_se", "use_use_RR_mean", "use_use_RR_se", 
                                  "use_mean_mean", "use_mean_se", "use_use_mean", "use_use_se")
for (g in 1:length(select_number)){
  for (h in 1:length(si)){
    for (i in 1:length(h2s)){
      for (j in 1:length(effective_marker_sizes)){
        A = results_cor_use[results_cor_use$number_offspring==select_number[g] & 
                              results_cor_use$si == si[h] & 
                              results_cor_use$h2s == h2s[i] & 
                              results_cor_use$effective_marker_sizes == effective_marker_sizes[j], ]
        results_cor_use_se[(g-1)*2*3*5 + (h-1)*3*5 + (i-1)*5 + j, ] = 
          c(select_number[g], si[h], h2s[i], effective_marker_sizes[j], 
            mean(A$use_mean_RR_cor), sd(A$use_mean_RR_cor)/sqrt(20), 
            mean(A$use_use_RR_cor), sd(A$use_use_RR_cor)/sqrt(20), 
            mean(A$use_mean_cor), sd(A$use_mean_cor)/sqrt(20), 
            mean(A$use_use_cor), sd(A$use_use_cor)/sqrt(20))
      }
    }
  }
}
results_cor_use_se_2 <- results_cor_use_se[results_cor_use_se$si==si[2], ]



p3 <- ggplot(results_cor_use_se_2, aes(number_offspring)) + 
  geom_point(aes(y=use_mean_mean, colour="predicted from family mean")) + 
  geom_errorbar(aes(ymin=use_mean_mean-use_mean_se, ymax=use_mean_mean+use_mean_se, 
                    colour="predicted from family mean"), width=20) + 
  geom_line(aes(y=use_mean_mean, colour="predicted from family mean"), 
            linewidth=0.2) + 
  geom_point(aes(y=use_use_mean, colour="predicted from family usefulness")) + 
  geom_errorbar(aes(ymin=use_use_mean-use_use_se, ymax=use_use_mean+use_use_se, 
                    colour="predicted from family usefulness"), width=20) + 
  geom_line(aes(y=use_use_mean, colour="predicted from family usefulness"), 
            linewidth=0.2) + 
  facet_grid(h2s~effective_marker_sizes) +
  # facet_wrap(~h2s, labeller = as_labeller(lbs, label_parsed)) + 
  xlab("number of offspring used") + 
  ylab("prediction accuracy of usefulness") + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  guides(color=guide_legend(title="prediction method", nrow=1)) + 
  # ggtitle("accuracy of predicting family usefulness of BV, BayesC") + 
  theme_minimal_grid(font_size=10) +
  theme(legend.position="bottom") 
  
save_plot(paste("view_correlation/plots/", "test2.pdf", sep=""),
          plot_grid(p3),
          base_width=6.5, base_height=4.33)




























