setwd("~/family_variance/")



marker_list <- readRDS("create_marker_list/phased_marker_list.rds")
marker_info_cM_parent <- readRDS("create_marker_list/marker_info_cM_parent.rds")
ind_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")
chrom_info <- readRDS("create_marker_list/chrom_info.rds")
chrom_table <- readRDS("create_marker_list/chrom_table.rds")



phasing_error <- c(0.05, 0.1, 0.15, 0.025, 0.01, 0.005)

introduce_error <- function(error, length){
  error_number = rpois(1, error*length)
  error_location = runif(error_number, 0, length)
  error_location = error_location[order(error_location)]
  return(error_location)
}

get_error_location <- function(error_location, chrom_info){
  # error_location: locations of error in cM
  # chrom_info: locations of marker in cM (with length of number of markers)
  # output: vector of true and false (with length of number of markers)
  # true indicates a switch between two markers, 
  # false not
  position = sapply(error_location, FUN=function(x){
    which(x<chrom_info)[1]
  })
  output = rep(F, length=length(chrom_info))
  for (i in 1:length(position)){
    output[position[i]:length(output)] = !output[position[i]:length(output)]
  }
  return(output)
}

set.seed(1)

for (h in 1:length(phasing_error)){
  marker_list_error = marker_list
  error_number = c()
  
  for (i in 1:length(ind_names)){
    error_list = lapply(chrom_table[, 3], introduce_error, error=phasing_error[h])
    error_numb = sapply(error_list, FUN=function(x){length(x)})
    error_number = cbind(error_number, error_numb)
    
    error_list_2 = vector(mode="list", length=length(chrom_names))
    for (j in 1:length(chrom_names)){
      if (length(error_list[[j]])==0){
        next 
      }
      error_list_2[[j]] = get_error_location(error_list[[j]], chrom_info[[j]])
      
      marker_list_error[[j]][error_list_2[[j]], , i] = marker_list[[j]][error_list_2[[j]], 2:1, i]
    }
  }
  
  saveRDS(error_number, paste("introduce_error/error_number_", 
                              phasing_error[h], 
                              ".rds", sep=""))
  saveRDS(marker_list_error, paste("introduce_error/marker_list_error_",
                                   phasing_error[h],
                                   ".rds", sep=""))
  
}


error_number_names <- list.files("introduce_error/", "error_number")

error_number <- lapply(error_number_names, FUN=function(x){
  readRDS(paste("introduce_error/", x, sep=""))
})



summary(as.vector(error_number[[1]]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   1.000   1.109   2.000   8.000 
summary(as.vector(error_number[[2]]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   1.000   2.000   2.229   3.000  13.000 
summary(as.vector(error_number[[3]]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   4.000   5.000   5.579   7.000  19.000





















