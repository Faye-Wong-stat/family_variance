setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)
sf <- c(0.8, 0.9, 0.95, 0.98)
b <- qnorm(sf, 0, 1)
si <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)

file_names_BV <- list.files("simulate_phenotypes_crosses/", pattern="BV.rds")
file_names_BV_fammean <- list.files("simulate_phenotypes_crosses/", pattern="BV_fammean.rds")
file_names_BV_famvar <- list.files("simulate_phenotypes_crosses/", pattern="BV_famvar.rds")
file_names_predY_RR_fammean <- list.files("simulate_phenotypes_crosses/", 
                                          pattern="predY_RR_fammean.rds")
file_names_predY_RR_famvar <- list.files("simulate_phenotypes_crosses/", 
                                         pattern="predY_RR_famvar.rds")
file_names_predY_fammean <- list.files("simulate_phenotypes_crosses/", 
                                       pattern="predY_fammean.rds")
file_names_predY_famvar <- list.files("simulate_phenotypes_crosses/", 
                                      pattern="predY_famvar.rds")
file_names <- gsub("_BV_fammean.rds", "", file_names_BV_fammean)
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
    vector(mode="list", length=nrow(parents_names))})
  for (k in 1:20){
    for (h in 1:nrow(parents_names)){
      parents_Z_noneffective[[j]][[k]][[h]] = Z[parents_names[h, ], -effective_marker_indices[[j]][, k]]
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

temp <- lapply(file_names_BV_fammean, FUN=function(x){
  readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
})
BV_fammean <- vector(mode="list", length=length(h2s))
BV_fammean <- lapply(BV_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(BV_fammean)){
  for (j in 1:length(BV_fammean[[i]])){
    BV_fammean[[i]][[j]] = temp[[1]][[i]][[j]]
    for (h in 2:length(file_names_BV_fammean)){
      BV_fammean[[i]][[j]] = rbind(BV_fammean[[i]][[j]], temp[[h]][[i]][[j]])
    }
  }
}

temp <- lapply(file_names_BV_famvar, FUN=function(x){
  readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
})
BV_famvar <- vector(mode="list", length=length(h2s))
BV_famvar <- lapply(BV_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(BV_famvar)){
  for (j in 1:length(BV_famvar[[i]])){
    BV_famvar[[i]][[j]] = temp[[1]][[i]][[j]]
    for (h in 2:length(file_names_BV_famvar)){
      BV_famvar[[i]][[j]] = rbind(BV_famvar[[i]][[j]], temp[[h]][[i]][[j]])
    }
  }
}

temp <- lapply(file_names_predY_RR_fammean, FUN=function(x){
  readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
})
predY_RR_fammean <- vector(mode="list", length=length(h2s))
predY_RR_fammean <- lapply(predY_RR_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(predY_RR_fammean)){
  for (j in 1:length(predY_RR_fammean[[i]])){
    predY_RR_fammean[[i]][[j]] = temp[[1]][[i]][[j]]
    for (h in 2:length(file_names_predY_RR_fammean)){
      predY_RR_fammean[[i]][[j]] = rbind(predY_RR_fammean[[i]][[j]], temp[[h]][[i]][[j]])
    }
  }
}

temp <- lapply(file_names_predY_RR_famvar, FUN=function(x){
  readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
})
predY_RR_famvar <- vector(mode="list", length=length(h2s))
predY_RR_famvar <- lapply(predY_RR_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(predY_RR_famvar)){
  for (j in 1:length(predY_RR_famvar[[i]])){
    predY_RR_famvar[[i]][[j]] = temp[[1]][[i]][[j]]
    for (h in 2:length(file_names_predY_RR_famvar)){
      predY_RR_famvar[[i]][[j]] = rbind(predY_RR_famvar[[i]][[j]], temp[[h]][[i]][[j]])
    }
  }
}

temp <- lapply(file_names_predY_fammean, FUN=function(x){
  readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
})
predY_fammean <- vector(mode="list", length=length(h2s))
predY_fammean <- lapply(predY_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(predY_fammean)){
  for (j in 1:length(predY_fammean[[i]])){
    predY_fammean[[i]][[j]] = temp[[1]][[i]][[j]]
    for (h in 2:length(file_names_predY_fammean)){
      predY_fammean[[i]][[j]] = rbind(predY_fammean[[i]][[j]], temp[[h]][[i]][[j]])
    }
  }
}

temp <- lapply(file_names_predY_famvar, FUN=function(x){
  readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
})
predY_famvar <- vector(mode="list", length=length(h2s))
predY_famvar <- lapply(predY_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(predY_famvar)){
  for (j in 1:length(predY_famvar[[i]])){
    predY_famvar[[i]][[j]] = temp[[1]][[i]][[j]]
    for (h in 2:length(file_names_predY_famvar)){
      predY_famvar[[i]][[j]] = rbind(predY_famvar[[i]][[j]], temp[[h]][[i]][[j]])
    }
  }
}

# create usefulness
BV_famuse <- vector(mode="list", length=length(si))
BV_famuse <- lapply(BV_famuse, FUN=function(x){
  vector(mode="list", length=length(h2s))
})
for (h in 1:length(si)){
  BV_famuse[[h]] = lapply(BV_famuse[[h]], FUN=function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
}
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      BV_famuse[[h]][[i]][[j]] = BV_fammean[[i]][[j]] + si[h] * sqrt(BV_famvar[[i]][[j]])
    }
  }
}
predY_RR_famuse <- vector(mode="list", length=length(si))
predY_RR_famuse <- lapply(predY_RR_famuse, FUN=function(x){
  vector(mode="list", length=length(h2s))
})
for (h in 1:length(si)){
  predY_RR_famuse[[h]] = lapply(predY_RR_famuse[[h]], FUN=function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
}
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      predY_RR_famuse[[h]][[i]][[j]] = 
        predY_RR_fammean[[i]][[j]] + si[h] * sqrt(predY_RR_famvar[[i]][[j]])
    }
  }
}
predY_famuse <- vector(mode="list", length=length(si))
predY_famuse <- lapply(predY_famuse, FUN=function(x){
  vector(mode="list", length=length(h2s))
})
for (h in 1:length(si)){
  predY_famuse[[h]] = lapply(predY_famuse[[h]], FUN=function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
}
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      predY_famuse[[h]][[i]][[j]] = 
        predY_fammean[[i]][[j]] + si[h] * sqrt(predY_famvar[[i]][[j]])
    }
  }
}




# predict fam mean bv
# predict fam sd bv with two methods 
# predict fam use with/without sd 
# corre between mean and sd*si true bv with slope 1 



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
cor_fammean_df <- data.frame(h2s = rep(paste("h2 = ", h2s, sep=""), each=(6*20)), 
                                effective_marker_sizes = 
                                  as.character(rep(rep(effective_marker_sizes, each=20), 5)), 
                                cor = rep(NA, (6*5*20)))
cor_fammean_df$effective_marker_sizes <- 
  factor(cor_fammean_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    cor_fammean_df$cor[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] = 
      cor_fammean[[i]][[j]]
  }
}
pdf(paste("view_usefulness/plots/", "fammean_cor.pdf", sep=""))
ggplot(cor_fammean_df, aes(as.numeric(effective_marker_sizes), cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~h2s) + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
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
# pdf(paste("view_usefulness_real_data/plots/", "fammean_cor_2.pdf", sep=""), width=10)
p1 <- ggplot(cor_fammean_df[cor_fammean_df$h2s %in% c("h2 = 0.9", "h2 = 0.5", "h2 = 0.1"), ], 
       aes(as.numeric(effective_marker_sizes), cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~h2s) + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  # ggtitle("accuracy of predicting family mean of BV, BayesC") + 
  theme_minimal_grid(font_size=10)
# dev.off()
save_plot(paste("view_usefulness_real_data/plots/", "fammean_cor_2.pdf", sep=""),
          p1, 
          base_width=6.5, base_height=2.17)



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
# cor_famsd_RR_df <- data.frame(h2s = rep(h2s, each=(6*20)), 
#                               effective_marker_sizes = 
#                                 as.character(rep(rep(effective_marker_sizes, each=20), 5)), 
#                               cor = rep(NA, (6*5*20)))
# cor_famsd_RR_df$effective_marker_sizes <- 
#   factor(cor_famsd_RR_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# for (i in 1:length(h2s)){
#   for (j in 1:length(effective_marker_sizes)){
#     cor_famsd_RR_df$cor[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] = 
#       cor_famsd_RR[[i]][[j]]
#   }
# }
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
# cor_famsd_df <- data.frame(h2s = rep(h2s, each=(6*20)), 
#                            effective_marker_sizes = 
#                              as.character(rep(rep(effective_marker_sizes, each=20), 5)), 
#                            cor = rep(NA, (6*5*20)))
# cor_famsd_df$effective_marker_sizes <- 
#   factor(cor_famsd_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# for (i in 1:length(h2s)){
#   for (j in 1:length(effective_marker_sizes)){
#     cor_famsd_df$cor[((i-1)*6*20 + (j-1)*20 + 1) : ((i-1)*6*20 + (j-1)*20 + 20)] = 
#       cor_famsd[[i]][[j]]
#   }
# }
cor_famsd_df <- data.frame(h2s = rep(paste("h2 = ", h2s, sep=""), each=(6*20)),
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

# # pdf(paste("view_usefulness_real_data/plots/", "famsd_cor_2.pdf", sep=""), width=10)
# p2.1 <- ggplot(cor_famsd_df[cor_famsd_df$h2s == 0.1, ], 
#        aes(as.numeric(effective_marker_sizes), y=cor)) +
#   geom_point() +
#   geom_smooth(method="loess") + 
#   ylim(0, 1) + 
#   xlab("number of causal loci") +
#   ylab("accuracy") + 
#   scale_x_continuous(breaks=1:6, labels=NULL) +
#   theme_minimal_grid(font_size=10) + 
#   theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
#         plot.title = element_text(hjust = 0.5)) + 
#   ggtitle("h2 = 0.1")
#   # ggtitle("accuracy of predicting family SD of BV, BayesC")
# # dev.off()
# p2.2 <- ggplot(cor_famsd_df[cor_famsd_df$h2s == 0.5, ], 
#                aes(as.numeric(effective_marker_sizes), y=cor)) +
#   geom_point() +
#   geom_smooth(method="loess") +
#   ylim(0, 1) + 
#   xlab("number of causal loci") +
#   ylab("accuracy") + 
#   scale_x_continuous(breaks=1:6, labels=NULL) +
#   theme_minimal_grid(font_size=10) + 
#   theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
#         plot.title = element_text(hjust = 0.5)) + 
#   ggtitle("h2 = 0.5")
# p2.3 <- ggplot(cor_famsd_df[cor_famsd_df$h2s == 0.9, ], 
#                aes(as.numeric(effective_marker_sizes), y=cor)) +
#   geom_point() +
#   geom_smooth(method="loess") +
#   ylim(0, 1) +
#   xlab("number of causal loci") +
#   ylab("accuracy") + 
#   scale_x_continuous(breaks=1:6, labels=NULL) +
#   theme_minimal_grid(font_size=10) + 
#   theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
#         plot.title = element_text(hjust = 0.5)) + 
#   ggtitle("h2 = 0.9")
# prow2 <- plot_grid(p2.1, p2.2, p2.3,  align="v", nrow=1)
# # pdf(paste("view_usefulness_real_data/plots/", "famsd_cor_two_way_2.pdf", sep=""), width=10)
# # gg_color_hue <- function(n){
# #   hues = seq(15, 375, length=n + 1)
# #   hcl(h=hues, l=65, c=100)[1:n]
# # }
# # cols <- gg_color_hue(length(levels(score_summ$error_level)))
# p3.1 <- ggplot(cor_famsd_df[cor_famsd_df$h2s == 0.1 , ], 
#        aes(as.numeric(effective_marker_sizes))) +
#   geom_point(aes(y=cor, colour="prediction from model")) +
#   geom_point(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents")) +
#   geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
#   geom_smooth(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents"), 
#               method="loess") +
#   ylim(-0.25, 1) +
#   xlab("number of causal loci") +
#   ylab("accuracy/correlation") +
#   scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
#   guides(color=guide_legend(nrow=1)) + 
#   theme(legend.position="bottom") +
#   # ggtitle("accuracy of predicting family SD of BV, BayesC")
#   theme_minimal_grid(font_size=10) + 
#   theme(strip.background = element_blank(),
#         strip.text.x = element_blank())
# # dev.off()
# p3.2 <- ggplot(cor_famsd_df[cor_famsd_df$h2s == 0.5, ], 
#                aes(as.numeric(effective_marker_sizes))) +
#   geom_point(aes(y=cor, colour="prediction from model")) +
#   geom_point(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents")) +
#   geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
#   geom_smooth(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents"), 
#               method="loess") +
#   theme(legend.position="bottom") +
#   ylim(-0.25, 1) +
#   xlab("number of causal loci") +
#   ylab("accuracy/correlation") +
#   scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
#   guides(color=guide_legend(nrow=1)) + 
#   theme_minimal_grid(font_size=10) + 
#   theme(strip.background = element_blank(),
#         strip.text.x = element_blank())

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
p1 <- ggplot(cor_fammean_df[cor_fammean_df$h2s %in% c(0.9, 0.5, 0.1), ], 
             aes(as.numeric(effective_marker_sizes), cor)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  facet_wrap(~h2s) + 
  xlab("number of causal loci") + 
  ylab("accuracy") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
  # ggtitle("accuracy of predicting family mean of BV, BayesC") + 
  theme_minimal_grid(font_size=10)
# dev.off()
save_plot(paste("view_usefulness/plots/", "fammean_cor_2.pdf", sep=""),
          p1, 
          base_width=6.5, base_height=2.17)

p3 <- ggplot(cor_famsd_df[cor_famsd_df$h2s %in% c("h2 = 0.9", "h2 = 0.5", "h2 = 0.1"), ], 
               aes(as.numeric(effective_marker_sizes))) +
  geom_point(aes(y=cor, colour="prediction from model")) +
  geom_point(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents")) +
  geom_smooth(aes(y=cor, colour="prediction from model"), method="loess") +
  geom_smooth(aes(y=cor_het_mar, colour="correlation with # of hetero markers of parents"), 
              method="loess") +
  facet_wrap(~h2s) + 
  ylim(-0.25, 1) +
  xlab("number of causal loci") +
  ylab("accuracy/correlation") +
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  guides(color=guide_legend(nrow=1)) + 
  theme_minimal_grid(font_size=10) + 
  theme(strip.background = element_blank())
# prow3 <- plot_grid(p3.1 + theme(legend.position="none"), 
#                    p3.2 + theme(legend.position="none"), 
#                    p3.3 + theme(legend.position="none"),  align="v", nrow=1)
# legend <- get_legend(p3.1)
# prow <- plot_grid(p2.1, p2.2, p2.3,
#                   p3.1 + theme(legend.position="none"), 
#                   p3.2 + theme(legend.position="none"), 
#                   p3.3 + theme(legend.position="none"), align="vh", ncol=3, nrow=2)
save_plot(paste("view_usefulness/plots/", "famsd_cor_two_way_2.pdf", sep=""), 
          plot_grid(p3 + theme(legend.position="none")), 
          base_width=6.5, base_height=2.17)



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
cor_famuse_df <- data.frame(si = rep(paste("si = ", round(si, 2), sep=""), each=(5*6*20)), 
                               h2s = rep(rep(paste("h2 = ", h2s, sep=""), each=(6*20)), 4), 
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
pdf(paste("view_usefulness/plots/", "famuse_cor.pdf", sep=""), width=14, height=10)
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
ggplot(cor_famuse_df[cor_famuse_df$h2s %in% c(0.9, 0.5, 0.1) & 
                       cor_famuse_df$si %in% c(1.4, 2.42), ], 
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



p4 <- ggplot(cor_famuse_df[cor_famuse_df$h2s %in% c("h2 = 0.9", "h2 = 0.5", "h2 = 0.1") & 
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
  guides(color=guide_legend(nrow=1)) + 
  # ggtitle("accuracy of predicting family usefulness of BV, BayesC") + 
  theme_minimal_grid(font_size=10)
save_plot(paste("view_usefulness_real_data/plots/", "famuse_cor_2.pdf", sep=""), 
          plot_grid(p4 + theme(legend.position="none")), 
          base_width=6.5, base_height=4.33)




# corre between usefulness and mean
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




# # corre between mean and sd*si true bv with slope 1 
# fammean_famsd <- vector(mode="list", length=length(si))
# fammean_famsd <- lapply(fammean_famsd, FUN=function(x){
#   vector(mode="list", length=length(effective_marker_sizes))
# })
# for (h in 1:length(si)){
#   for (j in 1:length(effective_marker_sizes)){
#     fammean_famsd[[h]][[j]] = matrix(NA, nrow=520*20, ncol=2)
#     colnames(fammean_famsd[[h]][[j]]) = c("mean", "sd")
#     fammean_famsd[[h]][[j]][, 1] = BV_fammean[[1]][[j]]
#     fammean_famsd[[h]][[j]][, 2] = si[h] * sqrt(BV_famvar[[1]][[j]])
#   }
# }
# fammean_famsd_df <- data.frame(si = rep(si, each=(6*520*20)), 
#                                effective_marker_sizes = 
#                                  as.character(rep(rep(effective_marker_sizes, each=(520*20)), 3)), 
#                                mean = rep(NA, 3*6*520*20), 
#                                sd = rep(NA, 3*6*520*20))
# fammean_famsd_df$effective_marker_sizes <- 
#   factor(fammean_famsd_df$effective_marker_sizes, levels=effective_marker_sizes)
# for (h in 1:length(si)){
#   for (j in 1:length(effective_marker_sizes)){
#     fammean_famsd_df$mean[((h-1)*6*520*20 + (j-1)*520*20 + 1) : 
#                            ((h-1)*6*520*20 + (j-1)*520*20 + 520*20)] = 
#       fammean_famsd[[h]][[j]][, 1]
#     fammean_famsd_df$sd[((h-1)*6*520*20 + (j-1)*520*20 + 1) : 
#                             ((h-1)*6*520*20 + (j-1)*520*20 + 520*20)] = 
#       fammean_famsd[[h]][[j]][, 2]
#   }
# }
# pdf(paste("view_usefulness_real_data/plots/", "fammean_famsd_cor.pdf", sep=""), width=10)
# ggplot(fammean_famsd_df, aes(x=mean, y=sd)) + 
#   geom_point() + 
#   geom_smooth(method="loess") + 
#   geom_abline(intercept=0, slope=1) + 
#   facet_grid(si~effective_marker_sizes) +
#   theme(legend.position="bottom") +
#   xlab("family mean") + 
#   ylab("family si*sd") +
#   # scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
#   ggtitle("correlation between family mean and si*sd of true BV")
# dev.off()
# pdf(paste("view_usefulness_real_data/plots/", "fammean_famsd_cor2.pdf", sep=""))
# for (h in 1:length(si)){
#   for (j in 1:length(effective_marker_sizes)){
#     print(
#       ggplot(as.data.frame(fammean_famsd[[h]][[j]]), aes(x=mean, y=sd)) + 
#         geom_point() + 
#         geom_smooth(method="loess") + 
#         geom_abline(intercept=0, slope=1) + 
#         xlab("family mean") + 
#         ylab("family si*sd") + 
#         ggtitle(paste("correlation between family mean and si*sd of true BV, si =", 
#                       si[h], 
#                       ", \nnumber of causal loci =", 
#                       effective_marker_sizes[j]))
#     )
#   }
# }
# dev.off()

# corre between mean and si true bv with slope 1 
fammean_famsd <- vector(mode="list", length=length(effective_marker_sizes))
for (j in 1:length(effective_marker_sizes)){
  fammean_famsd[[j]] = matrix(nrow=520*20, ncol=2)
  fammean_famsd[[j]][, 1] = BV_fammean[[1]][[j]]
  fammean_famsd[[j]][, 2] = sqrt(BV_famvar[[1]][[j]])
}
fammean_famsd_df <- data.frame(effective_marker_sizes=rep(effective_marker_sizes, each=520*20), 
                               mean=rep(NA, 6*520*20), 
                               sd=rep(NA, 6*520*20))
fammean_famsd_df$effective_marker_sizes <- 
  factor(fammean_famsd_df$effective_marker_sizes, levels=effective_marker_sizes)
for (j in 1:length(effective_marker_sizes)){
  fammean_famsd_df[((j-1)*520*20 + 1) : ((j-1)*520*20 + 520*20), 2:3] = 
    fammean_famsd[[j]]
}
pdf(paste("view_usefulness/plots/", "fammean_famsd.pdf", sep=""))
ggplot(fammean_famsd_df, aes(x=mean, y=sd)) + 
  geom_point() + 
  geom_abline(intercept=0, slope=1) + 
  facet_wrap(~effective_marker_sizes) + 
  xlab("family mean") + 
  ylab("family sd") + 
  ggtitle("family sd against mean")
dev.off()
pdf(paste("view_usefulness/plots/", "fammean_famsd2.pdf", sep=""))
for (j in 1:length(effective_marker_sizes)){
  print(
    ggplot(fammean_famsd_df[fammean_famsd_df$effective_marker_sizes==effective_marker_sizes[j], ], 
           aes(x=mean, y=sd)) + 
      geom_point() + 
      geom_abline(intercept=0, slope=1) + 
      xlab("family mean") + 
      ylab("family sd") + 
      ggtitle(paste("family sd against mean, number of causal loci = ", effective_marker_sizes[j]))
  )
}
dev.off()








# plotting family mean and offspring bvs
set.seed(1)
family_index <- sample(1:length(BV), 50, replace=F)
family_index <- family_index[order(family_index)]

BV_sample <- BV[family_index]
temp <- lapply(file_names_BV_fammean, FUN=function(x){
  readRDS(paste("simulate_phenotypes_crosses/", x, sep=""))
})
for (i in 1:length(temp)){
  temp[[i]] = temp[[i]][[1]]
}
BV_fammean_sample <- temp[family_index]


BV_BV_fammean_sample_df <- data.frame(family=rep(family_index, each=6*20*200), 
                                      effective_marker_sizes=
                                        rep(rep(effective_marker_sizes, each=20*200), 50), 
                                      pheno=rep(rep(1:20, each=200), 50*6), 
                                      BV_fammean=rep(NA, 50*6*20*200), 
                                      BV=rep(NA, 50*6*20*200))
BV_BV_fammean_sample_df$effective_marker_sizes <- 
  factor(BV_BV_fammean_sample_df$effective_marker_sizes, levels=effective_marker_sizes)
for (i in 1:length(BV_sample)){
  for (j in 1:length(BV_sample[[i]])){
    for (k in 1:20){
      BV_BV_fammean_sample_df$BV_fammean[((i-1)*6*20*200 + (j-1)*20*200 + (k-1)*200 + 1) : 
                                           ((i-1)*6*20*200 + (j-1)*20*200 + (k-1)*200 + 200)] = 
        rep(BV_fammean_sample[[i]][[j]][k], 200)
      BV_BV_fammean_sample_df$BV[((i-1)*6*20*200 + (j-1)*20*200 + (k-1)*200 + 1) : 
                                   ((i-1)*6*20*200 + (j-1)*20*200 + (k-1)*200 + 200)] = 
        BV_sample[[i]][[j]][, k]
    }
  }
}

pdf(paste("view_usefulness/plots/", "fammean_BV.pdf", sep=""))
ggplot(BV_BV_fammean_sample_df, aes(x=BV_fammean, y=BV)) + 
  geom_point() + 
  facet_wrap(~effective_marker_sizes) 
dev.off()
pdf(paste("view_usefulness/plots/", "fammean_BV2.pdf", sep=""))
for(j in 1:length(effective_marker_sizes)){
  print(
    ggplot(BV_BV_fammean_sample_df[
      BV_BV_fammean_sample_df$effective_marker_sizes == effective_marker_sizes[j], 
    ], aes(x=BV_fammean, y=BV)) + 
      geom_point() 
  )
}
dev.off()





# # usefulness correlation using rrBLUP
# cor_use_RR <- vector(mode="list", length=length(h2s))
# cor_use_RR <- lapply(cor_use_RR, FUN=function(x){
#   vector(mode="list", length=length(effective_marker_sizes))
# })
# 
# for (i in 1:length(BV_fammean)){
#   for (j in 1:length(BV_fammean[[1]])){
#     for (h in 1:20){
#       cor_use_RR[[i]][[j]][h] = 
#         cor(BV_fammean[[i]][[j]][, h] + si[2]*sqrt(BV_famvar[[i]][[j]][, h]), 
#             predY_RR_fammean[[i]][[j]][, h] + si[2]*sqrt(predY_RR_famvar[[i]][[j]][, h]))
#     }
#   }
# }
# usefulness_RR_cor <- data.frame(h2s = rep(h2s, each=(6*20)), 
#                              effective_marker_sizes = 
#                                as.character(rep(rep(effective_marker_sizes, each=20), 5)), 
#                              cor = rep(NA, (6*5*20)))
# usefulness_RR_cor$effective_marker_sizes <- 
#   factor(usefulness_RR_cor$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# corr <- c()
# for (i in 1:length(h2s)){
#   for (j in 1:length(effective_marker_sizes)){
#     corr = c(corr, cor_use_RR[[i]][[j]])
#   }
# }
# usefulness_RR_cor$cor <- corr
# 
# pdf(paste("view_usefulness_real_data/plots/", "useful_RR_cor.pdf", sep=""))
# ggplot(usefulness_RR_cor, aes(as.numeric(effective_marker_sizes), cor)) + 
#   geom_point() + facet_wrap(~h2s) + geom_smooth(method="loess") + xlab("effective marker sizes") + 
#   scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes))
# dev.off()
# pdf(paste("view_usefulness_real_data/plots/", "useful_RR_cor_xaxis_h2s.pdf", sep=""))
# ggplot(usefulness_RR_cor, aes(h2s, cor)) + 
#   geom_point() + facet_wrap(~effective_marker_sizes) + geom_smooth(method="loess")
# dev.off()
# 
# 
# 
# # usefulness correlation using bayesC
# cor_use <- vector(mode="list", length=length(h2s))
# cor_use <- lapply(cor_use, FUN=function(x){
#   vector(mode="list", length=length(effective_marker_sizes))
# })
# 
# for (i in 1:length(BV_fammean)){
#   for (j in 1:length(BV_fammean[[1]])){
#     for (h in 1:20){
#       cor_use[[i]][[j]][h] = 
#         cor(BV_fammean[[i]][[j]][, h] + si[2]*sqrt(BV_famvar[[i]][[j]][, h]), 
#             predY_fammean[[i]][[j]][, h] + si[2]*sqrt(predY_famvar[[i]][[j]][, h]))
#     }
#   }
# }
# usefulness_cor <- data.frame(h2s = rep(h2s, each=(6*20)), 
#                              effective_marker_sizes = 
#                                as.character(rep(rep(effective_marker_sizes, each=20), 5)), 
#                              cor = rep(NA, (6*5*20)))
# usefulness_cor$effective_marker_sizes <- 
#   factor(usefulness_cor$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# corr <- c()
# for (i in 1:length(h2s)){
#   for (j in 1:length(effective_marker_sizes)){
#     corr = c(corr, cor_use[[i]][[j]])
#   }
# }
# usefulness_cor$cor <- corr
# 
# pdf(paste("view_usefulness_real_data/plots/", "useful_cor.pdf", sep=""))
# ggplot(usefulness_cor, aes(as.numeric(effective_marker_sizes), cor)) + 
#   geom_point() + facet_wrap(~h2s) + geom_smooth(method="loess") + xlab("effective marker sizes") + 
#   scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes))
# dev.off()
# pdf(paste("view_usefulness_real_data/plots/", "useful_cor_xaxis_h2s.pdf", sep=""))
# ggplot(usefulness_cor, aes(h2s, cor)) + 
#   geom_point() + facet_wrap(~effective_marker_sizes) + geom_smooth(method="loess")
# dev.off()










