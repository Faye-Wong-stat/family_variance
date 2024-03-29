# select best parents for crosses based on true BV

setwd("~/family_variance/")
library(ggplot2)



indiv_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")

effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)

remove_marker_indices <- readRDS("simulate_phenotypes/remove_marker_indices.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")
Z <- readRDS("simulate_phenotypes/Z.rds")
alphas <- readRDS("simulate_phenotypes/alphas.rds")
Zalphas <- readRDS("simulate_phenotypes_real_data/Zalphas.rds")
# Zalphas_ls_mtx <-  readRDS("simulate_phenotypes/Zalphas_ls_mtx.rds")



best_parents <- data.frame(number_of_QTL=NA, replic=NA)
for (k in 1:20){
  best_parents[, ncol(best_parents)+1] = NA
  colnames(best_parents)[ncol(best_parents)] = paste("parent_", k, sep="")
}
for (j in 1:length(effective_marker_sizes)){
  for (k in 1:20){
    best_parents = rbind(best_parents, 
                         c(effective_marker_sizes[j], k, 
                           rownames(Z)[order(Zalphas[[j]][[k]][, 1], decreasing=T)[1:20]]))
  }
}

best_parents <- best_parents[-1, ]

j=1
k=1
Zalphas[[j]][[k]][order(Zalphas[[j]][[k]][, 1], decreasing=T)[1:20], 1]
# 16C045P030 09C048P003 10C126P002 10C144P001 91C248P003 16C049P026 04C041P006 
# 3.241007   3.241007   3.241007   3.241007   3.241007   3.057363   3.057363 
# 04C041P007 08C013P605 08C132P003 09C085P002 11C011P003 11C089P001 12C018P001 
# 3.057363   3.057363   3.057363   3.057363   3.057363   3.057363   3.057363 
# 12C023P604 12C029P601 12C053P601 12C061P605 12C072P606 12C075P001 
# 3.057363   3.057363   3.057363   3.057363   3.057363   3.057363 
best_parents[1:3, 1:5]
# number_of_QTL replic       parent_1   parent_2    parent_3
# 2             4      1     16C045P030 09C048P003  10C126P002
# 3             4      2     05C107P002 08C132P608 05C107P0021
# 4             4      3 Everbearing185   Hayazaki    PI551401

write.table(best_parents, "select_best_parents/best_parents.txt", row.names=F, col.names=F)
# write.table(medium_parents, "select_best_parents/medium_parents.txt", row.names=F, col.names=F)



# best_parents_long <- data.frame(number_of_QTL=NA, replic=NA, name=NA, BV=NA)
# for (j in 1:length(effective_marker_sizes)){
#   for (k in 1:20){
#     best_parents_long = rbind(best_parents_long, 
#                               data.frame(number_of_QTL=rep(effective_marker_sizes[j], 20), 
#                                          replic=rep(k, 20), 
#                                          name=rownames(Z)[order(Zalphas[[j]][[k]][, 1], 
#                                                                 decreasing=T)[1:20]], 
#                                          BV=Zalphas[[j]][[k]][, 1][order(Zalphas[[j]][[k]][, 1], 
#                                                                          decreasing=T)[1:20]]))
#   }
# }
# best_parents_long <- best_parents_long[-1, ]
# 
# best_parents_agg <- aggregate(best_parents_long$BV, list(best_parents_long$number_of_QTL, 
#                                                          best_parents_long$replic), mean)
# colnames(best_parents_agg) <- c("number_of_QTL", "replic", "mean")
# best_parents_agg$sd <- aggregate(best_parents_long$BV, list(best_parents_long$number_of_QTL, 
#                                                          best_parents_long$replic), sd)[, 3]
# 
# medium_parents_long <- data.frame(number_of_QTL=NA, replic=NA, name=NA, BV=NA)
# for (j in 1:length(effective_marker_sizes)){
#   for (k in 1:20){
#     medium_parents_long = rbind(medium_parents_long, 
#                               data.frame(number_of_QTL=rep(effective_marker_sizes[j], 20), 
#                                          replic=rep(k, 20), 
#                                          name=rownames(Z)[order(Zalphas[[j]][[k]][, 1], 
#                                                                 decreasing=T)[259:278]], 
#                                          BV=Zalphas[[j]][[k]][, 1][order(Zalphas[[j]][[k]][, 1], 
#                                                                          decreasing=T)[259:278]]))
#   }
# }
# medium_parents_long <- medium_parents_long[-1, ]
# 
# medium_parents_agg <- aggregate(medium_parents_long$BV, list(medium_parents_long$number_of_QTL, 
#                                                          medium_parents_long$replic), mean)
# colnames(medium_parents_agg) <- c("number_of_QTL", "replic", "mean")
# medium_parents_agg$sd <- aggregate(medium_parents_long$BV, list(medium_parents_long$number_of_QTL, 
#                                                             medium_parents_long$replic), sd)[, 3]
# 
# all_parents_long <- data.frame(number_of_QTL=NA, replic=NA, name=NA, BV=NA)
# for (j in 1:length(effective_marker_sizes)){
#   for (k in 1:20){
#     all_parents_long = rbind(all_parents_long, 
#                                 data.frame(number_of_QTL=rep(effective_marker_sizes[j], 537), 
#                                            replic=rep(k, 537), 
#                                            name=rownames(Z), 
#                                            BV=Zalphas[[j]][[k]][, 1]))
#   }
# }
# all_parents_long <- all_parents_long[-1, ]
# 
# all_parents_agg <- aggregate(all_parents_long$BV, list(all_parents_long$number_of_QTL, 
#                                                              all_parents_long$replic), mean)
# colnames(all_parents_agg) <- c("number_of_QTL", "replic", "mean")
# all_parents_agg$sd <- aggregate(all_parents_long$BV, list(all_parents_long$number_of_QTL, 
#                                                                 all_parents_long$replic), sd)[, 3]
# 
# 
# 
# parents_agg <- aggregate(best_parents_agg$mean, list(best_parents_agg$number_of_QTL), var)
# colnames(parents_agg) <- c("number_of_QTL", "var_of_mean_of_BV")
# parents_agg$var_of_sd_of_BV <- aggregate(best_parents_agg$sd, 
#                                         list(best_parents_agg$number_of_QTL), var)[, 2]
# parents_agg$parents_type <- "best parents"
# 
# parents_agg_01 <- aggregate(medium_parents_agg$mean, list(medium_parents_agg$number_of_QTL), var)
# colnames(parents_agg_01) <- c("number_of_QTL", "var_of_mean_of_BV")
# parents_agg_01$var_of_sd_of_BV <- aggregate(medium_parents_agg$sd, 
#                                         list(medium_parents_agg$number_of_QTL), var)[, 2]
# parents_agg_01$parents_type <- "medium parents"
# 
# parents_agg <- rbind(parents_agg, parents_agg_01)
# 
# parents_agg_02 <- aggregate(all_parents_agg$mean, list(all_parents_agg$number_of_QTL), var)
# colnames(parents_agg_02) <- c("number_of_QTL", "var_of_mean_of_BV")
# parents_agg_02$var_of_sd_of_BV <- aggregate(all_parents_agg$sd, 
#                                            list(all_parents_agg$number_of_QTL), var)[, 2]
# parents_agg_02$parents_type <- "all parents"
# 
# parents_agg <- rbind(parents_agg, parents_agg_02)
# 
# parents_agg$ratio <- parents_agg$var_of_sd_of_BV / parents_agg$var_of_mean_of_BV
# 
# parents_agg <- parents_agg[order(parents_agg$number_of_QTL, parents_agg$parents_type), ]
# parents_agg$number_of_QTL <- factor(parents_agg$number_of_QTL, levels=unique(parents_agg$number_of_QTL))
# 
# saveRDS(parents_agg, "select_best_parents/parents_agg.rds")
# 
# 
# 
# pdf("select_best_parents/plots/ratio_var_sd_var_mean.pdf")
# ggplot(parents_agg, aes(parents_type, ratio)) + 
#   geom_point() + 
#   facet_wrap(~number_of_QTL) + 
#   xlab("parents type") +
#   ylab("ratio of var sd and var mean") + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   ggtitle("ratio between var sd and var mean")
# dev.off()
# 
# pdf("select_best_parents/plots/ratio_var_sd_var_mean_QTL64.pdf")
# ggplot(parents_agg[parents_agg$number_of_QTL=="64", ], aes(parents_type, ratio)) + 
#   geom_point() + 
#   xlab("parents type") +
#   ylab("ratio of var sd and var mean") + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   ggtitle("ratio between var sd and var mean")
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# dim(Z)[1]/2 + 10
# dim(Z)[1]/2 - 10
# # 278.5
# # 258.5
# 
# medium_parents <- data.frame(number_of_QTL=NA, replic=NA)
# for (k in 1:20){
#   medium_parents[, ncol(medium_parents)+1] = NA
#   colnames(medium_parents)[ncol(medium_parents)] = paste("parent_", k, sep="")
# }
# for (j in 1:length(effective_marker_sizes)){
#   for (k in 1:20){
#     medium_parents = rbind(medium_parents, 
#                          c(effective_marker_sizes[j], k, 
#                            rownames(Z)[order(Zalphas[[j]][[k]][, 1], decreasing=T)[259:278]]))
#   }
# }
# 
# medium_parents <- medium_parents[-1, ]
# 
# Zalphas[[j]][[k]][order(Zalphas[[j]][[k]][, 1], decreasing=T)[259:278], 1]
# #  16C100P023  16C097P031  16C096P034  16C100P015  16C020P038  16C037P035 
# # -0.05443627 -0.05530602 -0.11920996 -0.11942452 -0.13404536 -0.14831346 
# #  16C021P032  16C036P022  16C046P003  16C561P023  16C052P026  16C047P039 
# # -0.19343467 -0.19724721 -0.22584422 -0.26280881 -0.27877876 -0.28200641 
# #  16C092P037  16C030P007  16C050P005  16C064P017  16C015P034  16C048P002 
# # -0.31442645 -0.31582930 -0.31870986 -0.34927425 -0.35103144 -0.36852057 
# #  16C048P032  16C059P005 
# # -0.36852057 -0.38282837











