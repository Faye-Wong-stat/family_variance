setwd("~/family_variance/")
library(ggplot2)



effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)

file_names <- list.files("simulate_crosses_best_parents/offspring_family_info_true/")
info_df <- sapply(file_names, FUN=function(x){
  strsplit(x, "_")[[1]][1:4]
})
info_df <- t(info_df)
colnames(info_df) <- c("number_of_QTL", "replic", "parent1", "parent2")
info_df <- as.data.frame(info_df)
info_df[, 1:2] <- apply(info_df[, 1:2], 2, as.numeric)

info_df$BV_fammean <- NA
info_df$BV_famsd <- NA

for (i in 1:length(file_names)){
  if(paste(paste(info_df[i, 1:4], collapse="_"), "_.rds", sep="") != file_names[i]){
    print(c(i, file_names[i]))
    next
  }
  offspring = readRDS(paste("simulate_crosses_best_parents/offspring_family_info_true/", file_names[i], sep=""))
  info_df$BV_fammean[i] = offspring[[2]]
  info_df$BV_famsd[i] = sqrt(offspring[[3]])

}

# info_df <- info_df[!is.na(info_df$parent2), ]

info_df$number_of_QTL <- factor(info_df$number_of_QTL, levels=effective_marker_sizes)
info_df$BV_use <- info_df$BV_fammean + info_df$BV_famsd

saveRDS(info_df, "view_crosses_best_parents/info_df.rds")


