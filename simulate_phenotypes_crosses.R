# setwd("~/family_variance/")
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args))
}

offsprings_name <- args[1]
trait_number <- as.numeric(args[2])

# sessionInfo()
# 360-365





# file_names <- read.delim("simulate_crosses/file_names.txt", header=F)

# i <- 499
# offsprings_name <- file_names[i,]
temp <- readRDS(offsprings_name)

offsprings_name <- gsub(".rds", "", gsub("simulate_crosses/offspring_list_", "", offsprings_name))
# trait_number <- 1

# parents_names <- gsub(".rds", "", gsub("simulate_crosses/offspring_list_", "", file_names[, 1]))
# parents_names <- matrix(unlist(strsplit(parents_names, split="_")), ncol=2, byrow=T)
# 
# for (i in 1:nrow(parents_names)){
#   if (any( apply(parents_names[-i, ], 1, setequal, parents_names[i, ]) )){
#     print("duplicated parents")
#     break
#   }
# }


# load the parameters 
effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)

lmms <- readRDS("simulate_phenotypes/lmms.rds")
bayesC <- readRDS("simulate_phenotypes/bayesC.rds")
remove_marker_indices <- readRDS("simulate_phenotypes/remove_marker_indices.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")
alphas <- readRDS("simulate_phenotypes/alphas.rds")
# var_epsilon <- readRDS("codes/simulate_phenotypes2_var_epsilon.rds")




# load a family's offsprings genotypes 
offsprings_genotype <- c()
for (j in 2:29){
  offsprings_genotype = rbind(offsprings_genotype, temp[[j]])
}
temp <- t(offsprings_genotype)
offsprings_genotype <- matrix(0, nrow=nrow(temp), ncol=ncol(temp))
offsprings_genotype[temp=="0|0"] <- -1
offsprings_genotype[temp=="1|1"] <- 1
offsprings_genotype <- offsprings_genotype[, !remove_marker_indices]
dim(offsprings_genotype)
# [1]   200 30899

# calculate offspring true BV
offsprings_Z <- vector(mode="list", length=length(alphas))
for (i in 1:length(offsprings_Z)){
  offsprings_Z[[i]] = offsprings_genotype
}
dim(offsprings_Z[[i]])
# [1]   200 30899

offsprings_Zalpha <- vector(mode="list", length=length(alphas))
for (i in 1:length(offsprings_Zalpha)){
  offsprings_Zalpha[[i]] = offsprings_Z[[i]] %*% alphas[[i]][[trait_number]]
    # matrix(alphas[[i]][[trait_number]], nrow=ncol(offsprings_Z[[i]]), ncol=1)
}
dim(offsprings_Zalpha[[i]])
# [1] 200   1

offsprings_BV <- vector(mode="list", length=length(h2s))
for (i in 1:length(h2s)){
  offsprings_BV[[i]] = offsprings_Zalpha
}

# calculate the family mean 
offsprings_BV_fammean <- vector(mode="list", length=length(h2s))
offsprings_BV_fammean <- lapply(offsprings_BV_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_BV_fammean[[i]][[j]] = mean(offsprings_BV[[i]][[j]])
  }
}
# length(offsprings_BV_fammean[[i]][[j]])
# [1] 1
# calculate the family variance 
offsprings_BV_famvar <- vector(mode="list", length=length(h2s))
offsprings_BV_famvar <- lapply(offsprings_BV_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_BV_famvar[[i]][[j]] = var(offsprings_BV[[i]][[j]])
  }
}
# length(offsprings_BV_famvar[[i]][[j]])
# [1] 1



# predict offsprings BV
offsprings_Z2 <- vector(mode="list", length=length(effective_marker_indices))
for (i in 1:length(offsprings_Z2)) {
  offsprings_Z2[[i]] = offsprings_genotype[, -effective_marker_indices[[i]][, trait_number]]
}
dim(offsprings_Z2[[i]])
# [1]   200 29875

offsprings_predY_RR <- vector(mode="list", length=length(lmms))
offsprings_predY_RR <- lapply(offsprings_predY_RR, FUN=function(x){
  vector(mode="list", length=length(lmms[[1]]))})
for (i in 1:length(offsprings_predY_RR)){
  for (j in 1:length(offsprings_predY_RR[[i]])){
    offsprings_predY_RR[[i]][[j]] = offsprings_Z2[[j]] %*% lmms[[i]][[j]][[trait_number]]$u
      # matrix(lmms[[i]][[j]][[trait_number]]$u, nrow=ncol(offsprings_Z2[[j]]))
  }
}

offsprings_predY_RR_fammean <- vector(mode="list", length=length(h2s))
offsprings_predY_RR_fammean <- lapply(offsprings_predY_RR_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_predY_RR_fammean[[i]][[j]] = mean(offsprings_predY_RR[[i]][[j]])
  }
}
offsprings_predY_RR_famvar <- vector(mode="list", length=length(h2s))
offsprings_predY_RR_famvar <- lapply(offsprings_predY_RR_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_predY_RR_famvar[[i]][[j]] = var(offsprings_predY_RR[[i]][[j]])
  }
}


offsprings_predY <- vector(mode="list", length=length(lmms))
offsprings_predY <- lapply(offsprings_predY, FUN=function(x){
  vector(mode="list", length=length(lmms[[1]]))})
for (i in 1:length(offsprings_predY)){
  for (j in 1:length(offsprings_predY[[i]])){
    offsprings_predY[[i]][[j]] = offsprings_Z2[[j]] %*% bayesC[[i]][[j]][[trait_number]]$ETA[[1]]$b
      # matrix(bayesC[[i]][[j]][[trait_number]]$ETA[[1]]$b, nrow=ncol(offsprings_Z2[[j]]))
  }
}

offsprings_predY_fammean <- vector(mode="list", length=length(h2s))
offsprings_predY_fammean <- lapply(offsprings_predY_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_predY_fammean[[i]][[j]] = mean(offsprings_predY[[i]][[j]])
  }
}
offsprings_predY_famvar <- vector(mode="list", length=length(h2s))
offsprings_predY_famvar <- lapply(offsprings_predY_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_predY_famvar[[i]][[j]] = var(offsprings_predY[[i]][[j]])
  }
}



# correlations between predicted BV and simulated BV
correlations_RR <- vector(mode="list", length=length(h2s))
correlations_RR <- lapply(correlations_RR, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(offsprings_predY_RR)){
  for (j in 1:length(offsprings_predY_RR[[i]])){
    correlations_RR[[i]][[j]] = cor(offsprings_BV[[i]][[j]], offsprings_predY_RR[[i]][[j]])
  }
}
names(correlations_RR) <- h2s
for (i in 1:length(correlations_RR)){
  names(correlations_RR[[i]]) <- effective_marker_sizes
}

correlations_RR_df <- data.frame(h2s=as.character(rep(h2s, each=5)),
                              effective_marker_sizes=
                                as.character(rep(effective_marker_sizes, 3)),
                              cor=rep(NA, 15))
correlations_RR_df$effective_marker_sizes <-
  factor(correlations_RR_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))

corr <- c()
for (i in 1:length(correlations_RR)){
  for (j in 1:length(correlations_RR[[i]])){
    corr = c(corr, correlations_RR[[i]][[j]])
    # print(c(i, j))
  }
}
correlations_RR_df$cor <- corr

correlations <- vector(mode="list", length=length(h2s))
correlations <- lapply(correlations, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(offsprings_predY)){
  for (j in 1:length(offsprings_predY[[i]])){
    correlations[[i]][[j]] = cor(offsprings_BV[[i]][[j]], offsprings_predY[[i]][[j]])
  }
}
names(correlations) <- h2s
for (i in 1:length(correlations)){
  names(correlations[[i]]) <- effective_marker_sizes
}

correlations_df <- data.frame(h2s=as.character(rep(h2s, each=5)),
                              effective_marker_sizes=
                                as.character(rep(effective_marker_sizes, 3)),
                              cor=rep(NA, 15))
correlations_df$effective_marker_sizes <-
  factor(correlations_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))

corr <- c()
for (i in 1:length(correlations)){
  for (j in 1:length(correlations[[i]])){
    corr = c(corr, correlations[[i]][[j]])
    # print(c(i, j))
  }
}
correlations_df$cor <- corr



result_df <- data.frame(family=rep(offsprings_name, 15), 
                        h2s=as.character(rep(h2s, each=5)),
                        effective_marker_sizes=
                         as.character(rep(effective_marker_sizes, 3)),
                        trait_number=rep(trait_number, 15),
                        BV_mean=rep(NA, 15), 
                        BV_var=rep(NA, 15), 
                        predY_RR_mean=rep(NA, 15), 
                        predY_RR_var=rep(NA, 15), 
                        predY_mean=rep(NA, 15), 
                        predY_var=rep(NA, 15))
result_df$effective_marker_sizes <- 
  factor(result_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    result_df[(i-1)*5 + j, 5:10] = c(offsprings_BV_fammean[[i]][[j]], 
                                offsprings_BV_famvar[[i]][[j]], 
                                offsprings_predY_RR_fammean[[i]][[j]], 
                                offsprings_predY_RR_famvar[[i]][[j]], 
                                offsprings_predY_fammean[[i]][[j]], 
                                offsprings_predY_famvar[[i]][[j]])
  }
}



saveRDS(offsprings_BV, paste("simulate_phenotypes_crosses/", 
                             offsprings_name, 
                            "_BV", 
                            ".rds", sep=""))
# saveRDS(offsprings_BV_fammean, paste("simulate_phenotypes_crosses/", 
#                              offsprings_name, 
#                              "_BV_fammean", 
#                              ".rds", sep=""))
# saveRDS(offsprings_BV_famvar, paste("simulate_phenotypes_crosses/", 
#                              offsprings_name, 
#                              "_BV_famvar", 
#                              ".rds", sep=""))
saveRDS(offsprings_predY_RR, paste("simulate_phenotypes_crosses/", 
                                offsprings_name, 
                                "_predY_RR", 
                                ".rds", sep=""))
# saveRDS(offsprings_predY_RR_fammean, paste("simulate_phenotypes_crosses/", 
#                                         offsprings_name, 
#                                 "_predY_RR_fammean", 
#                                 ".rds", sep=""))
# saveRDS(offsprings_predY_RR_famvar, paste("simulate_phenotypes_crosses/", 
#                                         offsprings_name, 
#                                         "_predY_RR_famvar", 
#                                         ".rds", sep=""))
saveRDS(offsprings_predY, paste("simulate_phenotypes_crosses/",
                                offsprings_name,
                                "_predY",
                                ".rds", sep=""))
# saveRDS(offsprings_predY_fammean, paste("simulate_phenotypes_crosses/",
#                                         offsprings_name,
#                                         "_predY_fammean",
#                                         ".rds", sep=""))
# saveRDS(offsprings_predY_famvar, paste("simulate_phenotypes_crosses/",
#                                        offsprings_name,
#                                        "_predY_famvar",
#                                        ".rds", sep=""))

saveRDS(correlations_RR_df, paste("simulate_phenotypes_crosses/",
                            offsprings_name,
                            "_cor_RR_df",
                            ".rds", sep=""))
saveRDS(correlations_df, paste("simulate_phenotypes_crosses/",
                            offsprings_name,
                            "_cor_df",
                            ".rds", sep=""))
saveRDS(result_df, paste("simulate_phenotypes_crosses/",
                         offsprings_name,
                         "_result_df",
                         ".rds", sep=""))


pdf(paste("simulate_phenotypes_crosses/plots/",
          offsprings_name,
          "_cor_RR",
          ".pdf", sep=""))
ggplot(correlations_RR_df, aes(x=effective_marker_sizes, y=cor)) + geom_point() + facet_wrap(~h2s)
dev.off()
pdf(paste("simulate_phenotypes_crosses/plots/",
          offsprings_name,
          "_cor",
          ".pdf", sep=""))
ggplot(correlations_df, aes(x=effective_marker_sizes, y=cor)) + geom_point() + facet_wrap(~h2s)
dev.off()










# full model
# correlations_full <- vector(mode="list", length=length(h2s))
# correlations_full <- lapply(correlations_full, FUN=function(x){
#   vector(mode="list", length=length(effective_marker_sizes))})
# for (i in 1:length(offsprings_predY_full)){
#   for (j in 1:length(offsprings_predY_full[[i]])){
#     for (k in 1:ncol(offsprings_predY_full[[i]][[j]])){
#       correlations_full[[i]][[j]][k] =
#         cor(offsprings_BV[[i]][[j]][, k], offsprings_predY_full[[i]][[j]][, k])
#     }
#   }
# }
# names(correlations_full) <- h2s
# for (i in 1:length(correlations_full)){
#   names(correlations_full[[i]]) <- effective_marker_sizes
# }
# 
# correlations_full_df <- data.frame(h2s=as.character(rep(h2s, each=120)),
#                                    effective_marker_sizes=
#                                      as.character(rep(rep(effective_marker_sizes, each=20), 5)),
#                                    cor=rep(NA, 600))
# correlations_full_df$effective_marker_sizes <-
#   factor(correlations_full_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# 
# corr <- c()
# for (i in 1:length(correlations_full)){
#   for (j in 1:length(correlations_full[[i]])){
#     corr = c(corr, correlations_full[[i]][[j]])
#     # print(c(i, j))
#   }
# }
# correlations_full_df$cor <- corr

# correlations_fammean_full <- data.frame(h2s=as.character(rep(h2s, each=6)),
#                                         effective_marker_sizes=as.character(rep(effective_marker_sizes,
#                                                                                 5)),
#                                         cor=rep(NA, 30))
# correlations_fammean_full$effective_marker_sizes <-
#   factor(correlations_fammean_full$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# corr <- c()
# for (i in 1:length(h2s)){
#   for (j in 1:length(effective_marker_sizes)){
#     corr = c(corr, cor(offsprings_BV_fammean[[i]][[j]], offsprings_predY_fammean_full[[i]][[j]]))
#   }
# }
# correlations_fammean_full$cor <- corr
# correlations_famvar_full <- data.frame(h2s=as.character(rep(h2s, each=6)),
#                                        effective_marker_sizes=as.character(rep(effective_marker_sizes,
#                                                                                5)),
#                                        cor=rep(NA, 30))
# correlations_famvar_full$effective_marker_sizes <-
#   factor(correlations_famvar_full$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
# corr <- c()
# for (i in 1:length(h2s)){
#   for (j in 1:length(effective_marker_sizes)){
#     corr = c(corr, cor(offsprings_BV_famvar[[i]][[j]], offsprings_predY_fammean_full[[i]][[j]]))
#   }
# }
# correlations_famvar_full$cor <- corr



# pdf(paste("simulation_result/",
#           offsprings_name,
#           "_cor_full",
#           ".pdf", sep=""))
# ggplot(correlations_full_df, aes(x=effective_marker_sizes, y=cor)) + geom_point() + facet_wrap(~h2s)
# dev.off()
# pdf(paste("simulation_result/",
#           offsprings_name,
#           "_cor_fammean_full",
#           ".pdf", sep=""))
# ggplot(correlations_fammean_full, aes(x=effective_marker_sizes, y=cor)) +
#   geom_point() + facet_wrap(~h2s)
# dev.off()
# pdf(paste("simulation_result/",
#           offsprings_name,
#           "_cor_famvar_full",
#           ".pdf", sep=""))
# ggplot(correlations_famvar_full, aes(x=effective_marker_sizes, y=cor)) +
#   geom_point() + facet_wrap(~h2s)
# dev.off()










# offsprings_Y <- readRDS(paste("simulation_result/", 
#                               gsub("simulated_data/", "", gsub(".rds", "", offsprings_name)), 
#                               "_Y", 
#                               ".rds", sep=""))
# offsprings_predY <- readRDS(paste("simulation_result/", 
#                                   gsub("simulated_data/", "", gsub(".rds", "", offsprings_name)), 
#                                   "_predY", 
#                                   ".rds", sep=""))
# correlations <- readRDS(paste("simulation_result/", 
#                               gsub("simulated_data/", "", gsub(".rds", "", offsprings_name)), 
#                               "_cor", 
#                               ".rds", sep=""))



# try without removing the causal markers 
# have a positive control 
# replace model$fit with true marker effects 
# compare cor(GV, estimated_GV)
# should be 1

# try average the family first 











# for (i in 1:length(correlations)){
#   correlations_df[(i-1)*120, 1] = names(correlations)[i]
#   for (j in 1:length(correlations[[i]])){
#     correlations_df[(i-1)*120 + (j-1)*20, 2] = names(correlations[[i]])[j]
#     correlations_df[((i-1)*120 + (j-1)*20) : ((i-1)*120 + (j-1)*20), 2]
#     
#     correlations_df[((i-1)*120 + (j-1)*20 + 1) : ((i-1)*120 + (j-1)*20 + 20), ] = 
#       c(names(correlations)[i], names(correlations[[i]])[j])
#   }
# }




# offsprings <- vector(mode="list", length=4)
# names(offsprings) <- c("parents names", "genotype", "predicted Y", "Y")
# offsprings[[1]] <- temp[[1]]
# for (j in 2:29){
#   offsprings[[2]] = rbind(offsprings[[2]], temp[[j]])
# }
# temp <- t(offsprings[[2]])
# offsprings[[2]] <- matrix(0, nrow=nrow(temp), ncol=ncol(temp))
# offsprings[[2]][temp=="0|0"] <- -1
# offsprings[[2]][temp=="1|1"] <- 1
# 
# offsprings[[2]] <- offsprings[[2]][, !remove_marker_indices]
# 
# 
# 
# offsprings[[3]] <- vector(mode="list", length=30) 
# for (j in 1:length(h2s)){
#   for (i in 1:length(effective_marker_sizes)){
#     names(offsprings[[3]])[(j-1)*5 + i] = paste(h2s[j], effective_marker_sizes[i], sep="_")
#     offsprings[[3]][[(j-1)*5 + i]] = 
#   }
# }




