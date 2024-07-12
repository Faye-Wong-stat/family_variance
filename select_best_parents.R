# select best parents for crosses based on predicted BV

setwd("~/family_variance/")
# library(ggplot2)



indiv_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")

effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)

remove_marker_indices <- readRDS("simulate_phenotypes/remove_marker_indices.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")
Z <- readRDS("simulate_phenotypes/Z.rds")
alphas <- readRDS("simulate_phenotypes/alphas.rds")
Zalphas <- readRDS("simulate_phenotypes/Zalphas.rds")
bayesC <- readRDS("simulate_phenotypes/bayesC.rds")

# Zalphas_ls_mtx <-  readRDS("simulate_phenotypes/Zalphas_ls_mtx.rds")



best_parents <- data.frame(h2s=NA, effective_marker_sizes=NA, replic=NA)
for (k in 1:40){
  best_parents[, ncol(best_parents)+1] = NA
  colnames(best_parents)[ncol(best_parents)] = paste("parent_", k, sep="")
}
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      predY = Z[, -effective_marker_indices[[j]][, k]] %*% 
        bayesC[[i]][[j]][[k]]$ETA[[1]]$b
      best_parents = rbind(best_parents, 
                            c(h2s[i], effective_marker_sizes[j], k, 
                              rownames(predY)[order(predY[, 1], decreasing=T)[1:40]]))
    }
  }
}

best_parents <- best_parents[-1, ]

i=1
j=1
k=1
predY = Z[, -effective_marker_indices[[j]][, k]] %*% 
  bayesC[[i]][[j]][[k]]$ETA[[1]]$b
predY[rownames(predY)[order(predY[, 1], decreasing=T)[1:40]], ]
# 16C056P013 12C004P001 16C546P024 87C112P006 16C036P029 16C062P008 12C155P001 
# 3.906551   3.785600   3.570153   3.394738   3.339412   3.312591   3.235645 
# 11C158P001 16C036P012 16C056P012 12C186P001 07C132P003 05C205P002 16C056P005 
# 3.222633   3.139962   2.979139   2.960802   2.957838   2.952749   2.895504 
# 97C207P003 16C061P020 12C059P602 12C071P001 11C103P001 16C036P022 12C083P003 
# 2.878256   2.878155   2.876393   2.861147   2.851466   2.840743   2.835140 
# 03C114P003 16C059P013 16C036P017 16C060P007 98C152P008 16C061P023 12C001P614 
# 2.796261   2.780907   2.729451   2.712650   2.683945   2.666673   2.602325 
# 12C143P001 11C036P601 12C089P002 10C167P004 16C055P010 12C180P004 00C295P002 
# 2.592333   2.577038   2.564352   2.546493   2.527103   2.524685   2.513309 
# 07C148P003 16C020P011 05C109P002 03C163P001 16C036P016 
# 2.509325   2.488332   2.476135   2.450194   2.444524 
best_parents[1:3, ]
# h2s effective_marker_sizes replic   parent_1   parent_2   parent_3   parent_4
# 2 0.8                      4      1 16C056P013 12C004P001 16C546P024 87C112P006
# 3 0.8                      4      2 16C556P001 16C092P026 35C093P011 16C072P008
# 4 0.8                      4      3 16C507P040 16C100P031 16C536P027 96C042P601
# parent_5   parent_6   parent_7   parent_8   parent_9  parent_10  parent_11
# 2 16C036P029 16C062P008 12C155P001 11C158P001 16C036P012 16C056P012 12C186P001
# 3 16C536P001 01C206P005 16C052P025 16C052P026 01C138P001 16C556P014 16C014P006
# 4 16C033P008 16C033P038 16C086P011 16C536P016 16C033P013 88C066P616 11C092P001
# parent_12  parent_13  parent_14  parent_15  parent_16  parent_17  parent_18
# 2 07C132P003 05C205P002 16C056P005 97C207P003 16C061P020 12C059P602 12C071P001
# 3 16C536P016 16C020P009 16C039P016 16C052P005 82C014P603 16C052P019 16C049P033
# 4 11C061P001 16C102P033 16C104P026 16C086P019 06C217P001 16C029P018 16C033P031
# parent_19  parent_20  parent_21  parent_22  parent_23  parent_24  parent_25
# 2 11C103P001 16C036P022 12C083P003 03C114P003 16C059P013 16C036P017 16C060P007
# 3 16C092P004 16C020P016 16C069P008 16C556P037 16C536P032 16C536P008 16C536P027
# 4 11C147P001 16C033P016 11C089P001 16C015P026 16C104P022 12C031P003 88C047P002
# parent_26  parent_27  parent_28  parent_29  parent_30  parent_31  parent_32
# 2 98C152P008 16C061P023 12C001P614 12C143P001 11C036P601 12C089P002 10C167P004
# 3 16C039P040 16C556P038 16C019P014 88C066P610 16C052P023 16C052P018 16C021P039
# 4 88C070P613      MSU61 16C033P028 01C206P005  MarysPeak    FVC1158 16C090P033
# parent_33  parent_34  parent_35  parent_36  parent_37  parent_38  parent_39
# 2 16C055P010 12C180P004 00C295P002 07C148P003 16C020P011 05C109P002 03C163P001
# 3 89C259P032 16C102P033 11C195P003 16C092P037 02C132P006 70C027P103 16C089P035
# 4 16C102P003 77C030P605 ORUS140776 16C040P016 16C536P032 16C100P020 11C105P001
# parent_40
# 2 16C036P016
# 3 16C042P014
# 4 11C036P002

write.table(best_parents, "select_best_parents/best_parents.txt", row.names=F, col.names=F)
# best_parents <- read.delim("select_best_parents/best_parents.txt", sep=" ", header=F)
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











