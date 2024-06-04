# setwd("~/family_variance/")
# library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args))
}

offsprings_name <- args[1]



# file_names <- read.delim("simulate_gametes/file_names.txt", header=F)
# i <- 1
# offsprings_name <- file_names[i,]

temp <- readRDS(offsprings_name)

offsprings_name <- gsub(".rds", "", gsub("simulate_gametes/", "", offsprings_name))



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)

remove_marker_indices <- readRDS("simulate_phenotypes/remove_marker_indices.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")
alphas <- readRDS("simulate_phenotypes/alphas.rds")
lmms <- readRDS("simulate_phenotypes/lmms.rds")
bayesC <- readRDS("simulate_phenotypes/bayesC.rds")





# load a family's offsprings genotypes 
offsprings_genotype <- c()
for (j in 2:29){
  offsprings_genotype = rbind(offsprings_genotype, temp[[j]])
}
offsprings_genotype <- t(offsprings_genotype)
offsprings_genotype <- offsprings_genotype[, !remove_marker_indices]
dim(offsprings_genotype)
# [1]   200 30899



# calculate gametes true BV
offsprings_Z <- vector(mode="list", length=length(alphas))
for (i in 1:length(offsprings_Z)){
  offsprings_Z[[i]] = offsprings_genotype
}
dim(offsprings_Z[[i]])
# [1]   200 30899

offsprings_Zalpha <- vector(mode="list", length=length(alphas))
for (i in 1:length(alphas)){
  offsprings_Zalpha[[i]] = matrix(NA, nrow=nrow(offsprings_Z[[i]]), ncol=length(alphas[[i]]))
  for (j in 1:length(alphas[[i]])){
    offsprings_Zalpha[[i]][, j] = offsprings_Z[[i]] %*% alphas[[i]][[j]]
  }
}
dim(offsprings_Zalpha[[i]])
# [1] 200   20

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
    offsprings_BV_fammean[[i]][[j]] = apply(offsprings_BV[[i]][[j]], 2, mean)
  }
}
# length(offsprings_BV_fammean[[i]][[j]])
# [1] 20
# calculate the family variance 
offsprings_BV_famvar <- vector(mode="list", length=length(h2s))
offsprings_BV_famvar <- lapply(offsprings_BV_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_BV_famvar[[i]][[j]] = apply(offsprings_BV[[i]][[j]], 2, var)
  }
}
# length(offsprings_BV_famvar[[i]][[j]])
# [1] 20



# calculate gametes predicted BV
offsprings_Z2 <- vector(mode="list", length=length(effective_marker_indices))
offsprings_Z2 <- lapply(offsprings_Z2, function(x){
  vector(mode="list", length=ncol(effective_marker_indices[[1]]))
})
for (i in 1:length(offsprings_Z2)) {
  for (j in 1:length(offsprings_Z2[[i]]))
  offsprings_Z2[[i]][[j]] = offsprings_genotype[, -effective_marker_indices[[i]][, j]]
}
dim(offsprings_Z2[[i]][[j]])
# [1] 200 29875

offsprings_predY_RR <- vector(mode="list", length=length(lmms))
offsprings_predY_RR <- lapply(offsprings_predY_RR, FUN=function(x){
  vector(mode="list", length=length(lmms[[1]]))})
for (i in 1:length(offsprings_predY_RR)){
  for (j in 1:length(offsprings_predY_RR[[i]])){
    offsprings_predY_RR[[i]][[j]] = offsprings_Z2[[j]][[1]] %*% lmms[[i]][[j]][[1]]$u
    for (k in 2:length(lmms[[i]][[j]])){
      offsprings_predY_RR[[i]][[j]] = cbind(offsprings_predY_RR[[i]][[j]], 
                                            offsprings_Z2[[j]][[k]] %*% lmms[[i]][[j]][[k]]$u)
    }
  }
}

offsprings_predY <- vector(mode="list", length=length(bayesC))
offsprings_predY <- lapply(offsprings_predY, FUN=function(x){
  vector(mode="list", length=length(bayesC[[1]]))})
for (i in 1:length(offsprings_predY)){
  for (j in 1:length(offsprings_predY[[i]])){
    offsprings_predY[[i]][[j]] = offsprings_Z2[[j]][[1]] %*% bayesC[[i]][[j]][[1]]$ETA[[1]]$b
    for (k in 2:length(bayesC[[i]][[j]])){
      offsprings_predY[[i]][[j]] = cbind(offsprings_predY[[i]][[j]], 
                                         offsprings_Z2[[j]][[k]] %*% bayesC[[i]][[j]][[k]]$ETA[[1]]$b)
    }
  }
}

offsprings_predY_RR_fammean <- vector(mode="list", length=length(h2s))
offsprings_predY_RR_fammean <- lapply(offsprings_predY_RR_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_predY_RR_fammean[[i]][[j]] = apply(offsprings_predY_RR[[i]][[j]], 2, mean)
  }
}
offsprings_predY_RR_famvar <- vector(mode="list", length=length(h2s))
offsprings_predY_RR_famvar <- lapply(offsprings_predY_RR_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_predY_RR_famvar[[i]][[j]] = apply(offsprings_predY_RR[[i]][[j]], 2, var)
  }
}

offsprings_predY_fammean <- vector(mode="list", length=length(h2s))
offsprings_predY_fammean <- lapply(offsprings_predY_fammean, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_predY_fammean[[i]][[j]] = apply(offsprings_predY[[i]][[j]], 2, mean)
  }
}
offsprings_predY_famvar <- vector(mode="list", length=length(h2s))
offsprings_predY_famvar <- lapply(offsprings_predY_famvar, FUN=function(x){
  vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    offsprings_predY_famvar[[i]][[j]] = apply(offsprings_predY[[i]][[j]], 2, var)
  }
}



result_df <- data.frame(family=rep(offsprings_name, 3*5*20), 
                        h2s=as.character(rep(h2s, each=5*20)),
                        effective_marker_sizes=
                          as.character(rep(rep(effective_marker_sizes, each=20), 3)),
                        trait_number=rep(1:20, 3*5),
                        BV_mean=rep(NA, 3*5*20), 
                        BV_var=rep(NA, 3*5*20), 
                        predY_RR_mean=rep(NA, 3*5*20), 
                        predY_RR_var=rep(NA, 3*5*20), 
                        predY_mean=rep(NA, 3*5*20), 
                        predY_var=rep(NA, 3*5*20))
result_df$effective_marker_sizes <- 
  factor(result_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    result_df[((i-1)*5*20 + (j-1)*20 + 1) : ((i-1)*5*20 + (j-1)*20 + 20), 
              5:10] = cbind(offsprings_BV_fammean[[i]][[j]], 
                                     offsprings_BV_famvar[[i]][[j]], 
                                     offsprings_predY_RR_fammean[[i]][[j]], 
                                     offsprings_predY_RR_famvar[[i]][[j]], 
                                     offsprings_predY_fammean[[i]][[j]], 
                                     offsprings_predY_famvar[[i]][[j]])
  }
}



saveRDS(offsprings_BV, paste("simulate_phenotypes_gametes/", 
                             offsprings_name, 
                             "_BV", 
                             ".rds", sep=""))
saveRDS(offsprings_predY_RR, paste("simulate_phenotypes_gametes/", 
                                   offsprings_name, 
                                   "_predY_RR", 
                                   ".rds", sep=""))
saveRDS(offsprings_predY, paste("simulate_phenotypes_gametes/",
                                offsprings_name,
                                "_predY",
                                ".rds", sep=""))
saveRDS(result_df, paste("simulate_phenotypes_gametes/",
                         offsprings_name,
                         "_result_df",
                         ".rds", sep=""))
