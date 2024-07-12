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

# create nested list with the same architecture 
# for example, len=c(2, 3)
# length(list) -> 2
# length(list[[1]]) -> 3
create.list <- function(len){
  if(length(len) == 1){
    vector("list", len)
  } else {
    lapply(1:len[1], function(...) create.list(len[-1]))
  }
}
# credit: 
# https://stackoverflow.com/questions/17567172/nested-lists-how-to-define-the-size-before-entering-data



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s                    <- c(0.8, 0.5, 0.2)
sf                     <- c(0.8, 0.99)
b                      <- qnorm(sf, 0, 1)
si                     <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)
indiv_names            <- readRDS("create_marker_list/indiv_names.rds")

BV_mean_family    <- readRDS("extract_gamete_info/BV_mean_family.R")
BV_sd_family      <- readRDS("extract_gamete_info/BV_sd_family.R")
BV_use_family     <- readRDS("extract_gamete_info/BV_use_family.R")
predY_mean_family <- readRDS("extract_gamete_info/predY_mean_family.R")
predY_sd_family   <- readRDS("extract_gamete_info/predY_sd_family.R")
predY_use_family  <- readRDS("extract_gamete_info/predY_use_family.R")

best_parents <- read.delim("select_best_parents/best_parents.txt", sep=" ", header=F)
dim(best_parents)
# [1] 300  43



BV_mean_family2    <- create.list(c(3, 5, 20))
BV_sd_family2      <- create.list(c(3, 5, 20))
BV_use_family2     <- create.list(c(2, 3, 5, 20))
predY_mean_family2 <- create.list(c(3, 5, 20))
predY_sd_family2   <- create.list(c(3, 5, 20))
predY_use_family2  <- create.list(c(2, 3, 5, 20))

for (i in 1:3){
  for (j in 1:5){
    for (k in 1:20){
      best_parents_names = unlist(best_parents[best_parents$V1 == h2s[i] & 
                                          best_parents$V2 == effective_marker_sizes[j] & 
                                          best_parents$V3 == k, -c(1:3)])
      best_parents_indices = match(best_parents_names, indiv_names)
      
      BV_mean_family2[[i]][[j]][[k]]    = 
        BV_mean_family[[j]][[k]][best_parents_indices, best_parents_indices]
      BV_sd_family2[[i]][[j]][[k]]      = 
        BV_sd_family[[j]][[k]][best_parents_indices, best_parents_indices]
      predY_mean_family2[[i]][[j]][[k]] = 
        predY_mean_family[[i]][[j]][[k]][best_parents_indices, best_parents_indices]
      predY_sd_family2[[i]][[j]][[k]]   = 
        predY_sd_family[[i]][[j]][[k]][best_parents_indices, best_parents_indices]
      
      
      colnames(BV_mean_family2[[i]][[j]][[k]])    = best_parents_names
      rownames(BV_mean_family2[[i]][[j]][[k]])    = best_parents_names
      colnames(BV_sd_family2[[i]][[j]][[k]])      = best_parents_names
      rownames(BV_sd_family2[[i]][[j]][[k]])      = best_parents_names
      colnames(predY_mean_family2[[i]][[j]][[k]]) = best_parents_names
      rownames(predY_mean_family2[[i]][[j]][[k]]) = best_parents_names
      colnames(predY_sd_family2[[i]][[j]][[k]])   = best_parents_names
      rownames(predY_sd_family2[[i]][[j]][[k]])   = best_parents_names
      
      for (h in 1:2){
        BV_use_family2[[h]][[i]][[j]][[k]]    = 
          BV_use_family[[h]][[j]][[k]][best_parents_indices, best_parents_indices]
        predY_use_family2[[h]][[i]][[j]][[k]] = 
          predY_use_family[[h]][[i]][[j]][[k]][best_parents_indices, best_parents_indices]
        
        colnames(BV_use_family2[[h]][[i]][[j]][[k]])    = best_parents_names
        rownames(BV_use_family2[[h]][[i]][[j]][[k]])    = best_parents_names
        colnames(predY_use_family2[[h]][[i]][[j]][[k]]) = best_parents_names
        rownames(predY_use_family2[[h]][[i]][[j]][[k]]) = best_parents_names
      }
    }
  }
}





# 3*5*20
best_parent_mean <- data.frame(h2s = rep(h2s, each=5*20), 
                               effective_marker_sizes = rep(rep(effective_marker_sizes, each=20), 3), 
                               trait_number = rep(1:20, 3*5))
best_parent_mean <- cbind(best_parent_mean, matrix(NA, 300, 10))

# 2*3*5*20*10
best_pred_mean <- data.frame(si = rep(si, each=3*5*20*10), 
                             h2s = rep(rep(h2s, each=5*20*10), 2), 
                             effective_marker_sizes = rep(rep(effective_marker_sizes, each=20*10), 2*3), 
                             trait_number = rep(rep(1:20, each=10), 2*3*5), 
                             family      = NA, 
                             BV_mean     = NA, 
                             BV_sd       = NA, 
                             BV_use      = NA, 
                             predY_mean  = NA, 
                             predY_sd    = NA, 
                             predY_use   = NA)

for(i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      index = order(predY_mean_family[[i]][[j]][[k]], decreasing=T)
      nr_nc = sapply(index, rowcol, 1007)
       
      nr_nc_2 = as.matrix(nr_nc[, 1], ncol=1)
      l=2
      while (ncol(nr_nc_2) < 10){
        if (!any(nr_nc[, l] %in% nr_nc_2)){
          nr_nc_2 = cbind(nr_nc_2, nr_nc[, l])
        } 
        l = l+1
      }
      
      families = apply(nr_nc_2, 2, function(x){
        indiv_names[x]
      })
      families = t(families)
      best_parent_mean[(i-1)*5*20 + (j-1)*20 + k, 4:13] = 
        paste(families[, 1], families[, 2], sep="_")
      
      for (h in 1:length(si)){
        best_pred_mean[((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 1) : 
                         ((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 10), 5] = 
          paste(families[, 1], families[, 2], sep="_")
        best_pred_mean[((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 1) : 
                         ((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 10), 6:11] = 
          cbind(apply(nr_nc_2, 2, function(x){ BV_mean_family[[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ BV_sd_family[[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ BV_use_family[[h]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_mean_family[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_sd_family[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_use_family[[h]][[i]][[j]][[k]][x[1], x[2]] }))
      }
    }
  }
}



best_parent_use <- data.frame(si = rep(si, each=3*5*20), 
                              h2s = rep(rep(h2s, each=5*20), 2), 
                               effective_marker_sizes = rep(rep(effective_marker_sizes, each=20), 2*3), 
                               trait_number = rep(1:20, 2*3*5))
best_parent_use <- cbind(best_parent_use, matrix(NA, 600, 10))

# 2*3*5*20*10
best_pred_use <- data.frame(si = rep(si, each=3*5*20*10), 
                             h2s = rep(rep(h2s, each=5*20*10), 2), 
                             effective_marker_sizes = rep(rep(effective_marker_sizes, each=20*10), 2*3), 
                             trait_number = rep(rep(1:20, each=10), 2*3*5), 
                             family      = NA, 
                             BV_mean     = NA, 
                             BV_sd       = NA, 
                             BV_use      = NA, 
                             predY_mean  = NA, 
                             predY_sd    = NA,
                             predY_use   = NA)

for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        index = order(predY_use_family[[h]][[i]][[j]][[k]], decreasing=T)
        nr_nc = sapply(index, rowcol, 1007)
        
        nr_nc_2 = as.matrix(nr_nc[, 1], ncol=1)
        l=2
        while (ncol(nr_nc_2) < 10){
          if (!any(nr_nc[, l] %in% nr_nc_2)){
            nr_nc_2 = cbind(nr_nc_2, nr_nc[, l])
          } 
          l = l+1
        }
        
        families = apply(nr_nc_2, 2, function(x){
          indiv_names[x]
        })
        families = t(families)
        best_parent_use[(h-1)*3*5*20 + (i-1)*5*20 + (j-1)*20 + k, 5:14] = 
          paste(families[, 1], families[, 2], sep="_")
        
        best_pred_use[((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 1) : 
                         ((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 10), 5] = 
          paste(families[, 1], families[, 2], sep="_")
        best_pred_use[((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 1) : 
                        ((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 10), 6:11] = 
          cbind(apply(nr_nc_2, 2, function(x){ BV_mean_family[[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ BV_sd_family[[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ BV_use_family[[h]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_mean_family[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_sd_family[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_use_family[[h]][[i]][[j]][[k]][x[1], x[2]] }))
      }
    }
  }
}



best_parent_mean2 <- data.frame(h2s = rep(h2s, each=5*20), 
                               effective_marker_sizes = rep(rep(effective_marker_sizes, each=20), 3), 
                               trait_number = rep(1:20, 3*5))
best_parent_mean2 <- cbind(best_parent_mean2, matrix(NA, 300, 10))

# 2*3*5*20*10
best_pred_mean2 <- data.frame(si = rep(si, each=3*5*20*10), 
                             h2s = rep(rep(h2s, each=5*20*10), 2), 
                             effective_marker_sizes = rep(rep(effective_marker_sizes, each=20*10), 2*3), 
                             trait_number = rep(rep(1:20, each=10), 2*3*5), 
                             family      = NA, 
                             BV_mean     = NA, 
                             BV_sd       = NA,
                             BV_use      = NA, 
                             predY_mean  = NA, 
                             predY_sd    = NA,
                             predY_use   = NA)

for(i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      index = order(predY_mean_family2[[i]][[j]][[k]], decreasing=T)
      nr_nc = sapply(index, rowcol, 40)
      
      nr_nc_2 = as.matrix(nr_nc[, 1], ncol=1)
      l=2
      while (ncol(nr_nc_2) < 10){
        if (!any(nr_nc[, l] %in% nr_nc_2)){
          nr_nc_2 = cbind(nr_nc_2, nr_nc[, l])
        } 
        l = l+1
      }
      
      families = apply(nr_nc_2, 2, function(x){
        colnames(predY_mean_family2[[i]][[j]][[k]])[x]
      })
      families = t(families)
      best_parent_mean2[(i-1)*5*20 + (j-1)*20 + k, 4:13] = 
        paste(families[, 1], families[, 2], sep="_")
      
      for (h in 1:length(si)){
        best_pred_mean2[((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 1) : 
                         ((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 10), 5] = 
          paste(families[, 1], families[, 2], sep="_")
        best_pred_mean2[((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 1) : 
                          ((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 10), 6:11] = 
          cbind(apply(nr_nc_2, 2, function(x){ BV_mean_family2[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ BV_sd_family2[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ BV_use_family2[[h]][[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_mean_family2[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_sd_family2[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_use_family2[[h]][[i]][[j]][[k]][x[1], x[2]] }))
      }
    }
  }
}



best_parent_use2 <- data.frame(si = rep(si, each=3*5*20), 
                              h2s = rep(rep(h2s, each=5*20), 2), 
                              effective_marker_sizes = rep(rep(effective_marker_sizes, each=20), 2*3), 
                              trait_number = rep(1:20, 2*3*5))
best_parent_use2 <- cbind(best_parent_use2, matrix(NA, 600, 10))

# 2*3*5*20*10
best_pred_use2 <- data.frame(si = rep(si, each=3*5*20*10), 
                            h2s = rep(rep(h2s, each=5*20*10), 2), 
                            effective_marker_sizes = rep(rep(effective_marker_sizes, each=20*10), 2*3), 
                            trait_number = rep(rep(1:20, each=10), 2*3*5), 
                            family      = NA, 
                            BV_mean     = NA, 
                            BV_sd       = NA,
                            BV_use      = NA, 
                            predY_mean  = NA, 
                            predY_sd    = NA,
                            predY_use   = NA)

for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        index = order(predY_use_family2[[h]][[i]][[j]][[k]], decreasing=T)
        nr_nc = sapply(index, rowcol, 40)
        
        nr_nc_2 = as.matrix(nr_nc[, 1], ncol=1)
        l=2
        while (ncol(nr_nc_2) < 10){
          if (!any(nr_nc[, l] %in% nr_nc_2)){
            nr_nc_2 = cbind(nr_nc_2, nr_nc[, l])
          } 
          l = l+1
        }
        
        families = apply(nr_nc_2, 2, function(x){
          colnames(predY_use_family2[[h]][[i]][[j]][[k]])[x]
        })
        families = t(families)
        best_parent_use2[(h-1)*3*5*20 + (i-1)*5*20 + (j-1)*20 + k, 5:14] = 
          paste(families[, 1], families[, 2], sep="_")
        
        best_pred_use2[((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 1) : 
                        ((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 10), 5] = 
          paste(families[, 1], families[, 2], sep="_")
        best_pred_use2[((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 1) : 
                         ((h-1)*3*5*20*10 + (i-1)*5*20*10 + (j-1)*20*10 + (k-1)*10 + 10), 6:11] = 
          cbind(apply(nr_nc_2, 2, function(x){ BV_mean_family2[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ BV_sd_family2[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ BV_use_family2[[h]][[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_mean_family2[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_sd_family2[[i]][[j]][[k]][x[1], x[2]] }), 
                apply(nr_nc_2, 2, function(x){ predY_use_family2[[h]][[i]][[j]][[k]][x[1], x[2]] }))
      }
    }
  }
}



saveRDS(best_parent_mean, "select_best_parents_gametes/best_parent_mean.rds")
write.table(best_parent_mean, "select_best_parents_gametes/best_parent_mean.txt", 
            row.names=F, col.names=F)
saveRDS(best_parent_use, "select_best_parents_gametes/best_parent_use.rds")
write.table(best_parent_use, "select_best_parents_gametes/best_parent_use.txt", 
            row.names=F, col.names=F)

saveRDS(best_parent_mean2, "select_best_parents_gametes/best_parent_mean2.rds")
write.table(best_parent_mean2, "select_best_parents_gametes/best_parent_mean2.txt", 
            row.names=F, col.names=F)
saveRDS(best_parent_use2, "select_best_parents_gametes/best_parent_use2.rds")
write.table(best_parent_use2, "select_best_parents_gametes/best_parent_use2.txt", 
            row.names=F, col.names=F)

saveRDS(best_pred_mean, "select_best_parents_gametes/best_pred_mean.rds")
saveRDS(best_pred_use, "select_best_parents_gametes/best_pred_use.rds")
saveRDS(best_pred_mean2, "select_best_parents_gametes/best_pred_mean2.rds")
saveRDS(best_pred_use2, "select_best_parents_gametes/best_pred_use2.rds")


