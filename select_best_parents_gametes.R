# select_best_parents_gametes.R
setwd("~/family_variance/")




# extract the row and column index of an element 
# given the index of it within a square matrix
rowcol <- function(ix, n){
  # ix: the index of an element with the square matrix
  # n: number of rows and columns in the matrix
  nr = ix %% n 
  if (nr == 0){
    nr = n
    nc = ix %/% n 
  } else {
    nc = ix %/% n + 1
  }
  return(c(nr, nc))
}



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s                    <- c(0.8, 0.5, 0.2)
sf                     <- c(0.8, 0.99)
b                      <- qnorm(sf, 0, 1)
si                     <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)
indiv_names            <- readRDS("create_marker_list/indiv_names.rds")

predY_mean_family <- readRDS("extract_gamete_info/predY_mean_family.R")
predY_use_family <- readRDS("extract_gamete_info/predY_use_family.R")



best_parent_mean <- data.frame(h2s = rep(h2s, each=5*20), 
                               effective_marker_sizes = rep(rep(effective_marker_sizes, each=20), 3), 
                               trait_number = rep(1:20, 3*5))
best_parent_mean <- cbind(best_parent_mean, matrix(NA, 300, 500))

for(i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      index = order(predY_mean_family[[i]][[j]][[k]], decreasing=T)[1:500]
      nr_nc = sapply(index, rowcol, 1007)
      families = apply(nr_nc, 2, function(x){
        indiv_names[x]
      })
      families = t(families)
      best_parent_mean[(i-1)*5*20 + (j-1)*20 + k, 4:503] = paste(families[, 1], families[, 2], sep="_")
    }
  }
}

best_parent_use <- data.frame(si = rep(si, each=3*5*20), 
                              h2s = rep(rep(h2s, each=5*20), 2), 
                               effective_marker_sizes = rep(rep(effective_marker_sizes, each=20), 2*3), 
                               trait_number = rep(1:20, 2*3*5))
best_parent_use <- cbind(best_parent_use, matrix(NA, 600, 500))

for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        index = order(predY_use_family[[h]][[i]][[j]][[k]], decreasing=T)[1:500]
        nr_nc = sapply(index, rowcol, 1007)
        families = apply(nr_nc, 2, function(x){
          indiv_names[x]
        })
        families = t(families)
        best_parent_use[(h-1)*3*5*20 + (i-1)*5*20 + (j-1)*20 + k, 4:503] = 
          paste(families[, 1], families[, 2], sep="_")
      }
    }
  }
}

saveRDS(best_parent_mean, "select_best_parents_gametes/best_parent_mean.rds")
saveRDS(best_parent_use, "select_best_parents_gametes/best_parent_use.rds")
write.table(best_parent_mean, "select_best_parents_gametes/best_parent_mean.txt", row.names=F, col.names=F)
write.table(best_parent_use, "select_best_parents_gametes/best_parent_use.txt", row.names=F, col.names=F)

