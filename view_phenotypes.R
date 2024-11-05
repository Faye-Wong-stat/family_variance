setwd("~/family_variance/")
library(rrBLUP)
library(BGLR)
library(scales)
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 1024)

Z <- readRDS("simulate_phenotypes/Z.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")
alphas <- readRDS("simulate_phenotypes/alphas.rds")
var_Zalpha <- readRDS("simulate_phenotypes/var_Zalpha.rds")



# look at heterozygous QTL in the population
hetero <- data.frame(effective_marker_sizes = rep(effective_marker_sizes, each=20), 
                        trait_number = rep(1:20, 5))
hetero <- cbind(hetero, matrix(NA, nrow=5*20, ncol=nrow(Z)))
colnames(hetero)[3:ncol(hetero)] <- rownames(Z)
for (i in 1:length(effective_marker_indices)){
  for (j in 1:20){
    index = effective_marker_indices[[i]][, j]
    marker_matrix = Z[, index]
    hetero[(i-1)*20 + j, 3:ncol(hetero)] = apply(marker_matrix, 1, function(x){
      sum(x==0)
    })
  }
}

heterozygosity <- data.frame(effective_marker_sizes = rep(effective_marker_sizes, each=20), 
                     trait_number = rep(1:20, 5), 
                     number_hetero_QTL=NA)
for (i in 1:length(effective_marker_indices)){
  for (j in 1:20){
    heterozygosity$number_hetero_QTL[(i-1)*20 + j] = 
      apply(hetero[(i-1)*20 + j, 3:ncol(hetero)], 1, function(x){
        mean(x)
      })
  }
}
heterozygosity$effective_marker_sizes <- factor(heterozygosity$effective_marker_sizes, 
                                                levels=as.factor(effective_marker_sizes))



pdf("view_phenotypes/plots/heterozygosity.pdf", width=3.5, height=3.5)
ggplot(heterozygosity, aes(effective_marker_sizes, number_hetero_QTL)) + 
  # geom_point() + 
  geom_boxplot(outlier.size = 0.75, size=0.4) + 
  xlab("number of causal loci") + 
  ylab("average number of heterozygous QTL per genotype") + 
  expand_limits(y=0.1) + 
  # scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes))  + 
  scale_y_continuous(trans="log10", breaks=c(0.1, 0.25, 0.5, 1, seq(0, 10, 2), 
                                             seq(0, 100, 20), seq(0, 1000, 200))) + 
  theme_minimal_grid(font_size=7) 
dev.off()

pdf("view_phenotypes/plots/heterozygosity2.pdf")
for (i in 1:length(effective_marker_sizes)){
  print(
    hist(unlist(hetero[hetero$effective_marker_sizes==effective_marker_sizes[i], 3:ncol(hetero)]), 
         main=paste("number of heterozygous QTL per genotype\n", 
                    "number of causal loci: ", effective_marker_sizes[i]), 
         xlab="number of heterozygous QTL per genotype", 
         ylab="count")
  )
}
dev.off()



# look at allele frequency
alle_freq <- colMeans(Z+1)/2
minor_alle_freq <- ifelse(alle_freq<0.5, alle_freq, (1-alle_freq))
pdf("view_phenotypes/plots/allele_frequency.pdf")
hist(colMeans(Z+1)/2, 
     main="allele frequency",
     xlab="allele frequency", 
     ylab="count", 
     xlim=c(0,1))
dev.off()
pdf("view_phenotypes/plots/minor_allele_frequency.pdf")
hist(minor_alle_freq, 
     main="allele frequency",
     xlab="allele frequency", 
     ylab="count", 
     xlim=c(0,1))
dev.off()



# look at marker effect sizes over certain threshold 
sd_Zalpha <- sqrt(var_Zalpha)
sd_Zalpha_20 <- sd_Zalpha*0.2
sd_Zalpha_10 <- sd_Zalpha*0.1
sd_Zalpha_05 <- sd_Zalpha*0.05
sd_Zalpha_all <- list(sd_Zalpha_20, sd_Zalpha_10, sd_Zalpha_05)

effect_size <- vector(mode="list", length(effective_marker_sizes))
effect_size <- lapply(effect_size, function(x){vector(mode="list", length=20)})
for (i in 1:length(effective_marker_sizes)){
  effect_size[[i]] = lapply(effect_size[[i]], function(x){vector(mode="list", length=3)})
  for (j in 1:20){
    for (k in 1:3){
      effect_size[[i]][[j]][[k]] = 
        which(alphas[[i]][[j]]!=0 & abs(alphas[[i]][[j]])>sd_Zalpha_all[[k]][j, i])
    }
  }
}

heterozygosity2 <- data.frame(effective_marker_sizes = rep(effective_marker_sizes, each=20*3), 
                             trait_number = rep(rep(1:20, each=3), 5), 
                             effect_size_percent = rep(c(20, 10, 5), ), 
                             number_hetero_QTL=NA)
for (i in 1:length(effective_marker_indices)){
  for (j in 1:20){
    for (k in 1:3){
      index = effect_size[[i]][[j]][[k]]
      marker_matrix = Z[, index, drop=F]
      heterozygosity2$number_hetero_QTL[(i-1)*20*3 + (j-1)*3 + k] = 
        mean(apply(marker_matrix, 1, function(x){sum(x==0)}))
    }
  }
}
heterozygosity2$effective_marker_sizes <- factor(heterozygosity2$effective_marker_sizes, 
                                                levels=as.factor(effective_marker_sizes))
heterozygosity2$effect_size_percent <- factor(heterozygosity2$effect_size_percent, 
                                              levels=as.factor(c(20, 10, 5)))
# heterozygosity2$number_hetero_QTL <- heterozygosity2$number_hetero_QTL+1

pdf("view_phenotypes/plots/effect_size.pdf", width=3.5, height=3.5)
ggplot(heterozygosity2, aes(effective_marker_sizes, number_hetero_QTL, color=effect_size_percent)) + 
  # geom_point() + 
  geom_boxplot(position="dodge", outlier.size = 0.75, size=0.4) +
  xlab("number of causal loci") + 
  ylab("average number of large effect heterozygous QTL per genotype") + 
  expand_limits(y=0.1) + 
  scale_y_continuous(trans=pseudo_log_trans(base=10), breaks=c(0.25, 0.5, 1, seq(0, 10, 2), 
                                             seq(0, 100, 20), seq(0, 1000, 200))) + 
  guides(color=guide_legend(title="effect size over certain\npercent of BV SD")) + 
  guides(color=guide_legend(title="prediction method")) + 
  # ggtitle("accuracy of predicting family usefulness of BV, BayesC") + 
  theme_minimal_grid(font_size=7) +
  theme(legend.position="bottom") 
dev.off()

heterozygosity2_mean <- aggregate(heterozygosity2$number_hetero_QTL, 
          by=list(heterozygosity2$effective_marker_sizes, heterozygosity2$effect_size_percent), 
          FUN=mean)
colnames(heterozygosity2_mean) <- c("effective_marker_sizes", 
                                    "effect_size_percent", 
                                    "number_hetero_QTL")
heterozygosity2_mean
# effective_marker_sizes effect_size_percent number_hetero_QTL
# 1                       4                  20        1.13396226
# 2                      16                  20        3.42919563
# 3                      64                  20        7.57820258
# 4                     256                  20        5.56886792
# 5                    1024                  20        0.06355511
# 6                       4                  10        1.28018868
# 7                      16                  10        4.56713009
# 8                      64                  10       13.97845084
# 9                     256                  10       29.78758689
# 10                   1024                  10       20.34880834
# 11                      4                   5        1.36281033
# 12                     16                   5        5.01812314
# 13                     64                   5       17.89682224
# 14                    256                   5       54.75129096
# 15                   1024                   5      116.96499503










