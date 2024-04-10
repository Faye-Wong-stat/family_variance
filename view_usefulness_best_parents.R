setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)
sf <- c(0.8, 0.9, 0.95, 0.98)
b <- qnorm(sf, 0, 1)
si <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)
si <- round(si, 2)

file_names_BV <- list.files("simulate_crosses_best_parents/offspring_family_info_true/")
file_names_pred <- list.files("simulate_crosses_best_parents/offspring_family_info_pred/")



cor_best_par <- matrix(NA, nrow=5*6*20*190, ncol=13)
cor_best_par <- as.data.frame(cor_best_par)
colnames(cor_best_par) <- c("h2", "number_of_QTL", "replic", "family",
                            "BV_mean", "BV_var", "BV_sd",
                            "pred_RR_mean", "pred_RR_var", "pred_RR_sd",
                            "pred_mean", "pred_var", "pred_sd")
cor_best_par$h2 <- rep(h2s, 6*20*190)

for (i in 1:length(file_names_BV)){
  name_BV = strsplit(file_names_BV[i], "_")[[1]]
  name_pred = strsplit(file_names_pred[i], "_")[[1]]
  if (any(name_BV != name_pred)){
    print("not the same family")
    print(i)
    break
  }

  object_BV = readRDS(paste("simulate_crosses_best_parents/offspring_family_info_true/",
                         file_names_BV[i], sep=""))
  object_pred = readRDS(paste("simulate_crosses_best_parents/offspring_family_info_pred/",
                            file_names_BV[i], sep=""))

  cor_best_par[((i-1)*5+1) : (i*5), c("number_of_QTL", "replic", "family")] =
    t(matrix(c(name_BV[1:2], paste(name_BV[3:4], collapse="_")), nrow=3, ncol=5))
  cor_best_par[((i-1)*5+1) : (i*5), c("BV_mean")] = object_BV[[2]]
  cor_best_par[((i-1)*5+1) : (i*5), c("BV_var")] = object_BV[[3]]

  cor_best_par[((i-1)*5+1) : (i*5), c("pred_RR_mean", "pred_RR_var", "pred_mean", "pred_var")] =
    matrix(unlist(object_pred[c(2:3, 5:6)]), nrow=5, ncol=4)
}

cor_best_par[, "BV_sd"] = sqrt(cor_best_par[, "BV_var"])
cor_best_par[, "pred_RR_sd"] = sqrt(cor_best_par[, "pred_RR_var"])
cor_best_par[, "pred_sd"] = sqrt(cor_best_par[, "pred_var"])
cor_best_par$number_of_QTL <- as.numeric(cor_best_par$number_of_QTL)
cor_best_par$replic <- as.numeric(cor_best_par$replic)

n <- nrow(cor_best_par)
cor_best_par <- cor_best_par[rep(1:n, length(si)), ]
cor_best_par$si <- rep(si, each=n)

cor_best_par$BV_uc <- cor_best_par$BV_mean + cor_best_par$si*cor_best_par$BV_sd
cor_best_par$pred_RR_uc <- cor_best_par$pred_RR_mean + cor_best_par$si*cor_best_par$pred_RR_sd
cor_best_par$pred_uc <- cor_best_par$pred_mean + cor_best_par$si*cor_best_par$pred_sd

saveRDS(cor_best_par, "view_usefulness_best_parents/cor_best_par.rds")
# cor_best_par <- readRDS("view_usefulness_best_parents/cor_best_par.rds")




# cor_best_par_df <- data.frame(h2 = rep(h2s, each=(6*4*20)),
#                               number_of_QTL =
#                                 rep(rep(effective_marker_sizes, each=4*20), 5),
#                               si = rep(rep(si, each=20), 5*6), 
#                               replic = rep(1:20, 5*6*4), 
#                               cor_mean_RR = rep(NA, (6*5*4*20)),
#                               cor_mean = rep(NA, (6*5*4*20)), 
#                               cor_sd_RR = rep(NA, (6*4*5*20)), 
#                               cor_sd = rep(NA, (6*5*4*20)))
# for (i in 1:nrow(cor_best_par_df)){
#   cor_best_par_tb = cor_best_par[cor_best_par$h2 == cor_best_par_df$h2[i] & 
#                                  cor_best_par$number_of_QTL == cor_best_par_df$number_of_QTL[i] & 
#                                  cor_best_par$si == cor_best_par_df$si[i] & 
#                                  cor_best_par$replic == cor_best_par_df$replic[i], ]
#   cor_best_par_df[i, c("cor_mean_RR", "cor_mean", "cor_sd_RR", "cor_sd")] = 
#     c(cor(cor_best_par_tb$BV_mean, cor_best_par_tb$pred_RR_mean),
#       cor(cor_best_par_tb$BV_mean, cor_best_par_tb$pred_mean),
#       cor(cor_best_par_tb$BV_sd, cor_best_par_tb$pred_RR_sd),
#       cor(cor_best_par_tb$BV_sd, cor_best_par_tb$pred_sd))
#   print(i)
# }
# cor_best_par_df$number_of_QTL <- as.character(cor_best_par_df$number_of_QTL)
# cor_best_par_df$number_of_QTL <- 
#   factor(cor_best_par_df$number_of_QTL, levels = as.factor(effective_marker_sizes))
# 
# n <- nrow(cor_best_par_df)
# cor_best_par_df <- do.call("rbind", replicate(length(si), cor_best_par_df, F))
# cor_best_par_df$si <- rep(si, each=n)

cor_best_par_df <- data.frame(h2 = rep(h2s, each=(6*4*20)),
                              number_of_QTL =
                                rep(rep(effective_marker_sizes, each=4*20), 5),
                              si = rep(rep(si, each=20), 5*6),
                              replic = rep(1:20, 5*6*4),
                              cor_ucBV_ucpred = rep(NA, (6*5*4*20)),
                              cor_ucBV_meanpred = rep(NA, (6*5*4*20)),
                              cor_ucBV_ucRR = rep(NA, (6*4*5*20)),
                              cor_ucBV_meanRR = rep(NA, (6*5*4*20)))
for (i in 1:nrow(cor_best_par_df)){
  cor_best_par_tb = cor_best_par[cor_best_par$h2 == cor_best_par_df$h2[i] &
                                 cor_best_par$number_of_QTL == cor_best_par_df$number_of_QTL[i] &
                                 cor_best_par$si == cor_best_par_df$si[i] &
                                 cor_best_par$replic == cor_best_par_df$replic[i], ]
  cor_best_par_df[i, c("cor_ucBV_ucpred", "cor_ucBV_meanpred", "cor_ucBV_ucRR", "cor_ucBV_meanRR")] =
    c(cor(cor_best_par_tb$BV_uc, cor_best_par_tb$pred_uc),
      cor(cor_best_par_tb$BV_uc, cor_best_par_tb$pred_mean),
      cor(cor_best_par_tb$BV_uc, cor_best_par_tb$pred_RR_uc),
      cor(cor_best_par_tb$BV_uc, cor_best_par_tb$pred_RR_mean))
  # print(i)
}
cor_best_par_df$number_of_QTL <- factor(cor_best_par_df$number_of_QTL, 
                                        levels=effective_marker_sizes)
cor_best_par_df$h2 <- paste("h^2 == ", cor_best_par_df$h2, sep="")
cor_best_par_df$si <- paste("i == ", cor_best_par_df$si ,sep="")



p1 <- ggplot(cor_best_par_df[cor_best_par_df$h2 %in% c("h^2 == 0.9", "h^2 == 0.5", "h^2 == 0.1") & 
                               cor_best_par_df$si %in% c("i == 1.4", "i == 2.42"), ], 
               aes(x=as.numeric(number_of_QTL))) +
  geom_point(aes(y=cor_ucBV_meanpred, colour="predicted from family mean BV")) +
  geom_point(aes(y=cor_ucBV_ucpred, colour="predicted from family usefulness BV")) +
  geom_smooth(aes(y=cor_ucBV_meanpred, colour="predicted from family mean BV"), method="loess") + 
  geom_smooth(aes(y=cor_ucBV_ucpred, colour="predicted from family usefulness BV"), method="loess") + 
  facet_grid(si~h2, labeller = label_parsed) + 
  ylim(0, 1) + 
  xlab("number of causal loci") +
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  theme_minimal_grid(font_size=10) +
  theme(legend.position="bottom") 

save_plot(paste("view_usefulness_best_parents/plots/", "fam_cor_best_parents.pdf", sep=""), 
          plot_grid(p1), 
          base_width=6.5, base_height=4.33)






