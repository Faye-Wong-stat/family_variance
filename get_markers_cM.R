setwd("family_variance/")
library(ggplot2)



# load the data
marker_info <- read.csv("generate_vcffiles/marker_info.csv", row.names=1)
marker_matrix_valid <- read.csv("generate_vcffiles/marker_matrix_valid.csv", row.names=1)
marker_info_valid <- read.csv("generate_vcffiles/marker_info_valid.csv", row.names=1)
marker_info_parent <- read.csv("generate_vcffiles/marker_info_parent.csv", row.names=1)
genetic_map <- read.csv("provided_data/20A823_Genetic_Map.csv", row.names=1)

genetic_map$Marker <- gsub("-", ".", genetic_map$Marker)

# # mark sure marker_info, marker_info_valid, and marker_info_parent are all the same
# for (i in 1:ncol(marker_info)){
#     if(!identical(marker_info[, i], marker_info_valid[, i])){
#       print(i)
#     }
# }
# for (i in 1:ncol(marker_info)){
#   if(!identical(marker_info[, i], marker_info_parent[, i])){
#     print(i)
#   }
# }


# create marker_distance and marker_distance_list
# create lms
marker_distance <- genetic_map[(genetic_map$Marker %in% rownames(marker_matrix_valid)) & 
                                 genetic_map$Status=="mapped", 
                               c("Marker", "Group_Nom_Alt", "cM_adj", "Cam_Pos")]
rownames(marker_distance) <- marker_distance$Marker
# marker_distance <- marker_distance[, -1]
colnames(marker_distance) <- c("ID", "chrom", "cM", "pos")
marker_distance <- marker_distance[order(marker_distance$chrom, marker_distance$cM), ]
marker_distance <- marker_distance[complete.cases(marker_distance), ]

chrom_names <- unique(marker_distance$chrom)
marker_distance_list <- list()
for (i in 1:length(chrom_names)){
  marker_distance_list[[i]] = marker_distance[marker_distance$chrom==chrom_names[i], ]
}
names(marker_distance_list) <- chrom_names

lm_list <- list()
for (i in 1:length(chrom_names)){
  lm_list[[i]] = lm(cM ~ pos, data=marker_distance_list[[i]])
}
names(lm_list) <- chrom_names

pdf("get_markers_cM/plots/regress cM based on 1.2k mapped markers.pdf")
for (i in 1:length(chrom_names)){
  print(
    ggplot(data=marker_distance_list[[i]], aes(x=pos, y=cM)) + 
      geom_point() + 
      geom_smooth(method="lm", formula=y~x) + 
      ggtitle(paste("cM against position in chromosome", chrom_names[i], "before data tidying")) + 
      xlab("position")
  )
  # print(plot(marker_distance_list[[i]]$pos, marker_distance_list[[i]]$cM,
  #            main=paste("cM against position in chromosome", chrom_names[i]),
  #            xlab="position", ylab="cM"),
  #       abline(lm(cM ~ pos, data=marker_distance_list[[i]]), col="red"))
}
dev.off()



# based on the plot, invert some markers and remove outlier markers 
# create marker_distance_cleaned and marker_distance_list_cleaned
marker_distance_list_cleaned <- marker_distance_list

{
# inversion on chrom 1B
marker_distance_list_cleaned[["1B"]][marker_distance_list_cleaned[["1B"]]$pos < 5e+06, "cM"] <- 
  rev(marker_distance_list_cleaned[["1B"]][marker_distance_list_cleaned[["1B"]]$pos < 5e+06, "cM"])

# outlier on chrom 1D
marker_distance_list_cleaned[["1D"]] <- 
  marker_distance_list_cleaned[["1D"]][!(marker_distance_list_cleaned[["1D"]]$pos < 5e+06 & 
                                       marker_distance_list_cleaned[["1D"]]$cM < 60), ]

# outlier on chrom 2A
marker_distance_list_cleaned[["2A"]] <- 
  marker_distance_list_cleaned[["2A"]][!(marker_distance_list_cleaned[["2A"]]$pos > 1.8e+07 & 
                                       marker_distance_list_cleaned[["2A"]]$cM > 1350), ]
marker_distance_list_cleaned[["2A"]] <- 
  marker_distance_list_cleaned[["2A"]][!(marker_distance_list_cleaned[["2A"]]$pos > 2.4e+07), ]

# inversion and outlier on chrom 2B
marker_distance_list_cleaned[["2B"]] <- 
  marker_distance_list_cleaned[["2B"]][!(marker_distance_list_cleaned[["2B"]]$pos < 3e+06 & 
                                           marker_distance_list_cleaned[["2B"]]$cM > 1850), ]
marker_distance_list_cleaned[["2B"]][marker_distance_list_cleaned[["2B"]]$pos > 1.8e+07 & 
                                       marker_distance_list_cleaned[["2B"]]$pos < 2.2e+07, "cM"] <- 
  rev(marker_distance_list_cleaned[["2B"]][marker_distance_list_cleaned[["2B"]]$pos > 1.8e+07 & 
                                             marker_distance_list_cleaned[["2B"]]$pos < 2.2e+07, "cM"])

# inversion on chrom 2C
marker_distance_list_cleaned[["2C"]][marker_distance_list_cleaned[["2C"]]$pos < 1.8e+07, "cM"] <- 
  rev(marker_distance_list_cleaned[["2C"]][marker_distance_list_cleaned[["2C"]]$pos < 1.8e+07, "cM"])

# outlier on chrom 2D
marker_distance_list_cleaned[["2D"]] <- 
  marker_distance_list_cleaned[["2D"]][!(marker_distance_list_cleaned[["2D"]]$pos < 5.0e+06 & 
                                       marker_distance_list_cleaned[["2D"]]$cM < 1500), ]
marker_distance_list_cleaned[["2D"]] <- 
  marker_distance_list_cleaned[["2D"]][!(marker_distance_list_cleaned[["2D"]]$pos > 1.5e+07 & 
                                         marker_distance_list_cleaned[["2D"]]$cM > 1500), ]
marker_distance_list_cleaned[["2D"]] <- 
  marker_distance_list_cleaned[["2D"]][!(marker_distance_list_cleaned[["2D"]]$pos < 2e+06), ]

# outlier on chrom 3A
marker_distance_list_cleaned[["3A"]] <- 
  marker_distance_list_cleaned[["3A"]][!(marker_distance_list_cleaned[["3A"]]$pos > 1.5e+07), ]

# inversion and outliers on chrom 3B 
marker_distance_list_cleaned[["3B"]] <- 
  marker_distance_list_cleaned[["3B"]][!(marker_distance_list_cleaned[["3B"]]$pos > 2.5e+07 & 
                                           marker_distance_list_cleaned[["3B"]]$cM < 2300), ]
marker_distance_list_cleaned[["3B"]][marker_distance_list_cleaned[["3B"]]$pos < 2.5e+06, "cM"] <- 
  rev(marker_distance_list_cleaned[["3B"]][marker_distance_list_cleaned[["3B"]]$pos < 2.5e+06, "cM"])

# outliers on chrom 3C
marker_distance_list_cleaned[["3C"]] <- 
  marker_distance_list_cleaned[["3C"]][!(marker_distance_list_cleaned[["3C"]]$cM < 2550 & 
                                       marker_distance_list_cleaned[["3C"]]$pos > 2.5e+06), ]

# outliers on chrom 3D
marker_distance_list_cleaned[["3D"]] <- 
  marker_distance_list_cleaned[["3D"]][!(marker_distance_list_cleaned[["3D"]]$pos < 5e+06 & 
                                       marker_distance_list_cleaned[["3D"]]$cM < 2050), ]

# outliers on chrom 4A
marker_distance_list_cleaned[["4A"]] <- 
  marker_distance_list_cleaned[["4A"]][!(marker_distance_list_cleaned[["4A"]]$pos < 2.1e+07 & 
                                       marker_distance_list_cleaned[["4A"]]$cM < 3410), ]

# outliers on chrom 4B
marker_distance_list_cleaned[["4B"]] <- 
  marker_distance_list_cleaned[["4B"]][!(marker_distance_list_cleaned[["4B"]]$pos > 1.3e+07 & 
                                       marker_distance_list_cleaned[["4B"]]$cM > 3700), ]

# outliers on chrom 4D
marker_distance_list_cleaned[["4D"]] <- 
  marker_distance_list_cleaned[["4D"]][!(marker_distance_list_cleaned[["4D"]]$pos < 8e+06 & 
                                       marker_distance_list_cleaned[["4D"]]$cM > 3025), ]

# outliers on chrom 5A 
marker_distance_list_cleaned[["5A"]] <- 
  marker_distance_list_cleaned[["5A"]][!(marker_distance_list_cleaned[["5A"]]$cM > 4065), ]
marker_distance_list_cleaned[["5A"]] <- 
  marker_distance_list_cleaned[["5A"]][!(marker_distance_list_cleaned[["5A"]]$pos > 2.5e+07 & 
                                       marker_distance_list_cleaned[["5A"]]$cM < 3900), ]

# outliers on chrom 5B
marker_distance_list_cleaned[["5B"]] <- 
  marker_distance_list_cleaned[["5B"]][!(marker_distance_list_cleaned[["5B"]]$pos > 2.5e+07), ]

# outliers on chrom 5D
marker_distance_list_cleaned[["5D"]] <- 
  marker_distance_list_cleaned[["5D"]][!(marker_distance_list_cleaned[["5D"]]$pos < 1e+07 & 
                                       marker_distance_list_cleaned[["5D"]]$cM > 4250), ]

# outliers on chrom 6B
marker_distance_list_cleaned[["6B"]] <- 
  marker_distance_list_cleaned[["6B"]][!(marker_distance_list_cleaned[["6B"]]$pos > 2e+07 & 
                                       marker_distance_list_cleaned[["6B"]]$cM > 5450), ]
marker_distance_list_cleaned[["6B"]] <- 
  marker_distance_list_cleaned[["6B"]][!(marker_distance_list_cleaned[["6B"]]$pos > 4e+07), ]

# inversions on chrom 6C
marker_distance_list_cleaned[["6C"]] <- 
  marker_distance_list_cleaned[["6C"]][!(marker_distance_list_cleaned[["6C"]]$pos > 3e+07 & 
                                           marker_distance_list_cleaned[["6C"]]$cM < 4950), ]
marker_distance_list_cleaned[["6C"]] <- 
  marker_distance_list_cleaned[["6C"]][!(marker_distance_list_cleaned[["6C"]]$pos < 5e+06 & 
                                         marker_distance_list_cleaned[["6C"]]$cM < 4900), ]
marker_distance_list_cleaned[["6C"]][marker_distance_list_cleaned[["6C"]]$pos < 2.7e+07, "cM"] <- 
  rev(marker_distance_list_cleaned[["6C"]][marker_distance_list_cleaned[["6C"]]$pos < 2.7e+07, "cM"])
marker_distance_list_cleaned[["6C"]][marker_distance_list_cleaned[["6C"]]$pos > 2.7e+07, "cM"] <- 
  rev(marker_distance_list_cleaned[["6C"]][marker_distance_list_cleaned[["6C"]]$pos > 2.7e+07, "cM"])

# outliers on chrom 6D
marker_distance_list_cleaned[["6D"]] <- 
  marker_distance_list_cleaned[["6D"]][!(marker_distance_list_cleaned[["6D"]]$cM > 5800 & 
                                       marker_distance_list_cleaned[["6D"]]$pos < 5e+06), ]

# inversion on chrom 7A
marker_distance_list_cleaned[["7A"]][marker_distance_list_cleaned[["7A"]]$pos > 1e+07 & 
                                       marker_distance_list_cleaned[["7A"]]$pos < 1.8e+7, "cM"] <- 
  rev(marker_distance_list_cleaned[["7A"]][marker_distance_list_cleaned[["7A"]]$pos > 1e+07 & 
                                             marker_distance_list_cleaned[["7A"]]$pos < 1.8e+7, "cM"])

# outliers on chrom 7B
marker_distance_list_cleaned[["7B"]] <- 
  marker_distance_list_cleaned[["7B"]][!(marker_distance_list_cleaned[["7B"]]$pos > 2e+07 & 
                                       marker_distance_list_cleaned[["7B"]]$cM > 6350), ]

# inversion on chrom 7C
marker_distance_list_cleaned[["7C"]][marker_distance_list_cleaned[["7C"]]$pos > 9e+06 & 
                                       marker_distance_list_cleaned[["7C"]]$pos < 1.7e+07, "cM"] <- 
  rev(marker_distance_list_cleaned[["7C"]][marker_distance_list_cleaned[["7C"]]$pos > 9e+06 & 
                                             marker_distance_list_cleaned[["7C"]]$pos < 1.7e+07, "cM"])

# outliers on chrom 7D
marker_distance_list_cleaned[["7D"]] <- 
  marker_distance_list_cleaned[["7D"]][!(marker_distance_list_cleaned[["7D"]]$pos < 5e+06 & 
                                       marker_distance_list_cleaned[["7D"]]$cM < 6500), ]
marker_distance_list_cleaned[["7D"]] <- 
  marker_distance_list_cleaned[["7D"]][!(marker_distance_list_cleaned[["7D"]]$pos > 2.4e+07), ]
}

pdf("get_markers_cM/plots/regress loess cM cleaned.pdf")
for (i in 1:length(chrom_names)){
  print(
    ggplot(data=marker_distance_list_cleaned[[i]], aes(x=pos, y=cM)) + 
      geom_point() + 
      geom_smooth(method="loess") + 
      ggtitle(paste("cM against position in chromosome", chrom_names[i], "after data tidying")) + 
      xlab("position")
  )
}
dev.off()
pdf("get_markers_cM/plots/regress lm cM cleaned.pdf")
for (i in 1:length(chrom_names)){
  print(
    ggplot(data=marker_distance_list_cleaned[[i]], aes(x=pos, y=cM)) + 
      geom_point() + 
      geom_smooth(method="lm", formula=y~x) + 
      ggtitle(paste("cM against position in chromosome", chrom_names[i], "after data tidying")) + 
      xlab("position")
  )
}
dev.off()

marker_distance_cleaned <- marker_distance_list_cleaned[[1]]
for (i in 2:length(chrom_names)){
  marker_distance_cleaned = rbind(marker_distance_cleaned, marker_distance_list_cleaned[[i]])
}
marker_info_clean <- marker_info[marker_info$ID %in% rownames(marker_matrix_valid),
                                 c("ID", "CHROM", "POS")]

# plot the distribution of markers from two data sets 
pdf("get_markers_cM/plots/distributions of markers.pdf")
for (i in 1:length(chrom_names)){
  print(
    ggplot(data=marker_distance_list_cleaned[[i]], aes(x=pos, fill="mapped markers")) + 
      geom_histogram(bins=50, position="identity", alpha=0.5) +
      geom_histogram(data=marker_info_clean[marker_info_clean$Chromosome==chrom_names[i], ], 
                     mapping=aes(x=POS, fill="phasing markers"), 
                     bins=50, position="identity", alpha=0.5) +
      ggtitle(paste("distribution of marker positions on chromosome", chrom_names[i])) +
      scale_color_manual(name="type",
                         breaks=c("mapped markers", "phasing markers"),
                         values=c("mapped markers"="black", "phasing markers"="red"))
  )
}
dev.off()
# 3A, 5B, 6B don't have good coverage



# create loesses
loess_list <- list()
for (i in 1:length(chrom_names)){
  if (chrom_names[i] %in% c("3A", "5B", "6B")){
    loess_list[[i]] = lm(cM ~ pos, data=marker_distance_list_cleaned[[i]])
  } else if(chrom_names[i]=="6C"){
    loess_list[[i]] = loess(cM ~ pos, data=marker_distance_list_cleaned[[i]], span=0.35)
  } else{
    loess_list[[i]] = loess(cM ~ pos, data=marker_distance_list_cleaned[[i]])
  }
}
names(loess_list) <- chrom_names
for (i in 1:length(chrom_names)){
  if (chrom_names[i] %in% c("3A", "5B", "6B")){
    marker_distance_list_cleaned[[i]]$loess_pred = loess_list[[i]]$fitted.values
  } else{
    marker_distance_list_cleaned[[i]]$loess_pred = loess_list[[i]]$fitted
  }
}

pdf("get_markers_cM/plots/regress loess cM cleaned 2.pdf")
for (i in 1:length(chrom_names)){
  print(
    ggplot(data=marker_distance_list_cleaned[[i]], aes(x=pos, y=cM, color="mapped markers")) + 
      geom_point() + 
      geom_line(aes(x=pos, y=loess_pred, color="fitted line")) + 
      scale_color_manual(name="type",
                         breaks=c("mapped markers", "fitted line"),
                         values=c("mapped markers"="black", "fitted line"="blue")) + 
      ggtitle(paste("cM against position in chromosome", chrom_names[i], "after data tidying")) + 
      xlab("position")
  )
}
dev.off()

{
# pdf("plots/get_markers_diff_dist/regress loess cM based on 1.2k mapped markers cleaned adjusted.pdf")
# for (i in 1:length(chrom_names)){
#   if (chrom_names[i]=="6C"){
#     print(
#       ggplot(data=marker_distance_list_cleaned[[i]], aes(x=pos, y=cM)) + 
#         geom_point() + 
#         geom_smooth(method="loess", span=0.25) + 
#         ggtitle(paste("cM against position in chromosome", chrom_names[i], 
#         "after data tidying, span=0.25")) + 
#         xlab("position")
#     )
#   } else{
#     print(
#       ggplot(data=marker_distance_list_cleaned[[i]], aes(x=pos, y=cM)) + 
#         geom_point() + 
#         geom_smooth(method="loess") + 
#         ggtitle(paste("cM against position in chromosome", chrom_names[i], "after data tidying")) + 
#         xlab("position")
#     )
#   }
# }
# dev.off()
}



# predict cM based on loess and lm
marker_info_cM_parent <- marker_info_parent
marker_info_cM_valid <- marker_info_valid
marker_info_cM_parent$cM <- NA
marker_info_cM_valid$cM <- NA

for (i in 1:length(chrom_names)){
  marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[i], "cM"] = 
    predict(loess_list[[i]], 
            newdata=data.frame(pos=marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[i], "POS"]))
}
sum(is.na(marker_info_cM_parent$cM))
# [1] 2371

pdf("get_markers_cM/plots/predict cM based on 1.2k mapped markers.pdf")
for (i in 1:length(chrom_names)){
  print(
    ggplot(data=marker_distance_list_cleaned[[i]], aes(x=pos, y=cM, color="training markers")) +
      geom_point(alpha=0.75) +
      geom_line(aes(x=pos, y=loess_pred, color="fitted line")) +
      geom_point(data=marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[i], ],
                 aes(x=POS, y=cM, color="predicted markers"), 
                 alpha=0.05) +
      ggtitle(paste("cM against position in chromosome", chrom_names[i])) +
      scale_color_manual(name="type",
                         breaks=c("training markers", "fitted line", "predicted markers"),
                         values=c("training markers"="black", "fitted line"="blue", "predicted markers"="red"))
  )
}
dev.off()

marker_info_cM_valid$cM <- marker_info_cM_parent$cM

write.csv(marker_info_cM_parent, "get_markers_cM/marker_info_cM_parent.csv")
write.csv(marker_info_cM_valid, "get_markers_cM/marker_info_cM_valid.csv")



# marker_info_cM_parent <- t(marker_info_cM_parent)
# distance <- list()
# for (i in 1:length(chrom_names)){
#   distance[[i]] = dist(marker_info_cM_parent[10, marker_info_cM_parent["CHROM", ]==chrom_names[i]])
# }



# head(marker_info_cM_valid[rownames(marker_distance[marker_distance$chrom==chrom_names[1], ]), ])
# #              CHROM     POS ID REF ALT QUAL FILTER INFO FORMAT       cM
# # AX.184070632    1A 2,138,519 NA   G   N   99   PASS   NA     GT 48.22367
# # AX.184907259    1A 1967827 NA   G   N   99   PASS   NA     GT 46.75316
# # AX.184022785    1A 2827535 NA   A   N   99   PASS   NA     GT 54.15951
# # AX.184629331    1A 2890915 NA   T   N   99   PASS   NA     GT 54.70553
# # AX.184025896    1A 3298167 NA   T   N   99   PASS   NA     GT 58.21398
# # AX.184070649    1A 3845701 NA   C   N   99   PASS   NA     GT 62.93097
# marker_info[marker_info$probe_id %in% c("AX.184070632", "AX.184907259", "AX.184022785"), ]
# #           probe_id class  chrom   pos.x hit_count hit_score ref_hit ref_site ref_nt allele_status chrom_status Chromosome   pos.y
# # 4221  AX.184022785   PHR Fvb1-4 2827535         1      5459  Fvb1-4  3330715      A  allele_match  chrom_match         1A 3330715
# # 6454  AX.184070632   PHR Fvb1-4 2138519         1      5332  Fvb1-4   957152      G  allele_match  chrom_match         1A  957152
# # 35380 AX.184907259   PHR Fvb1-4 1967827         1      5539  Fvb1-4  1106771      G  allele_match  chrom_match         1A 1106771
# head(marker_distance[marker_distance$chrom==chrom_names[1], ])
# #              chrom     cM      pos
# # AX.184070632    1A 725.81 80,992,877
# # AX.184907259    1A 726.82 81142496
# # AX.184022785    1A 748.72 83366440
# # AX.184629331    1A 748.72 83426473
# # AX.184025896    1A 753.28 83800080
# # AX.184070649    1A 759.93 84371943



# marker_info_clean <- marker_info[marker_info$probe_id %in% rownames(marker_info_cM_valid), 
#                                  c("probe_id", "Chromosome", "pos.x", "pos.y", "ref_site")]
# marker_info_clean2 <- marker_info_clean[marker_info_clean$probe_id %in% rownames(marker_distance), ]
# marker_info_clean2 <- cbind(marker_info_clean2, 
#                             marker_distance[match(marker_info_clean2$probe_id, rownames(marker_distance)), c("pos", "pos2")])
# colnames(marker_info_clean2)[6:7] <- c("Cam_Pos", "pos_adj")
# marker_info_clean2 <- marker_info_clean2[order(marker_info_clean2$Chromosome), ]
# pdf("plots/get_markers_diff_dist/three positions.pdf")
# for (i in 1:length(chrom_names)){
#   print(
#     plot(marker_info_clean2[marker_info_clean2$Chromosome==chrom_names[i], c(3:4, 6:7)], 
#          main=paste("three position columns in matrix info dataset for chromosome", chrom_names[i]))
#   )
# }
# dev.off()

