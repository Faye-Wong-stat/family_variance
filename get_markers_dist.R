setwd("~/family_variance/")
# library(vcfR)



distindex <- function(i, j, n){
  # i < j <= n
  # i:row, j:col, n:number of observations 
  return(n*(i-1) - i*(i-1)/2 + j - i)
}
rowcol <- function(ix, n){
  nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc=n-(2*n-nr+1)*nr/2+ix+nr
  return(c(nr,nc))
}



marker_info_cM_parent <- readRDS("get_markers_cM/marker_info_cM_parent.rds")
# marker_info_cM_valid <- read.csv("get_markers_cM/marker_info_cM_valid.csv")
# marker_matrix <- read.csv("cleaned_data/marker_matrix.csv", row.names=1, check.names=F)
# phased_marker_parent <- read.vcfR("phased_data/phased_parent.vcf")
# phased_marker_info_parent <- as.data.frame(phased_marker_parent@fix)
# phased_marker_matrix_parent <- as.data.frame(phased_marker_parent@gt)

# phased_marker_info_parent$FORMAT <- phased_marker_matrix_parent$FORMAT
# phased_marker_matrix_parent <- phased_marker_matrix_parent[, -1]
# rownames(phased_marker_info_parent) <- rownames(marker_matrix)
# rownames(phased_marker_matrix_parent) <- rownames(marker_matrix)
# for (i in 1:ncol(marker_info_cM_valid)){
#   print(identical(marker_info_cM_valid[, i], marker_info_cM_parent[, i]))
# }



marker_info_cM_parent2 <- marker_info_cM_parent[!is.na(marker_info_cM_parent$cM), ]
# marker_info_cM_valid <- marker_info_cM_valid[!is.na(marker_info_cM_valid$cM), ]
dim(marker_info_cM_parent)
# [1] 37282    10
dim(marker_info_cM_parent2)
# [1] 34911    10



chrom_names <- unique(marker_info_cM_parent[,"CHROM"])
cM_dist <- vector(mode="list", length=length(chrom_names))
names(cM_dist) <- chrom_names
for (j in 1:length(chrom_names)){
  cM_dist[[j]] = dist(marker_info_cM_parent2[marker_info_cM_parent2$CHROM==chrom_names[[j]], ]$cM)
}

for (j in 1:length(chrom_names)){
  print(max(cM_dist[[j]]))
}
# [1] 174.4135
# [1] 243.1175
# [1] 148.2739
# [1] 224.5724
# [1] 155.932
# [1] 280.6426
# [1] 300.3577
# [1] 240.2779
# [1] 251.6255
# [1] 274.5326
# [1] 213.508
# [1] 249.7291
# [1] 157.5015
# [1] 249.4875
# [1] 236.7943
# [1] 221.7733
# [1] 171.196
# [1] 218.7748
# [1] 200.9504
# [1] 234.787
# [1] 229.3853
# [1] 305.8853
# [1] 345.8841
# [1] 303.2631
# [1] 109.8548
# [1] 210.523
# [1] 105.9246
# [1] 195.8764



cM_dist_markers <- vector(mode="list", length=length(chrom_names))
names(cM_dist_markers) <- chrom_names
for (j in 1:length(chrom_names)){
  cM_Size = attr(cM_dist[[j]], "Size")
  
  dist_100 = matrix(NA, nrow=1, ncol=2)
  dist_10 = matrix(NA, nrow=1, ncol=2)
  dist_1 = matrix(NA, nrow=1, ncol=2)
  dist_0.1 = matrix(NA, nrow=1, ncol=2)
  dist_0.01 = matrix(NA, nrow=1, ncol=2)
  
  for (h in 1:length(cM_dist[[j]])){
    if (cM_dist[[j]][h] > 100){
      new_index = rowcol(h, cM_Size)
      if (!any(new_index %in% dist_100)){
        dist_100 = rbind(dist_100, new_index)
      }
    } else if(cM_dist[[j]][h] <= 100 & cM_dist[[j]][h] > 10){
      new_index = rowcol(h, cM_Size)
      if (!any(new_index %in% dist_10)){
        dist_10 = rbind(dist_10, new_index)
      }
    } else if(cM_dist[[j]][h] <= 10 & cM_dist[[j]][h] > 1){
      new_index = rowcol(h, cM_Size)
      if (!any(new_index %in% dist_1)){
        dist_1 = rbind(dist_1, new_index)
      }
    } else if(cM_dist[[j]][h] <= 1 & cM_dist[[j]][h] > 0.1){
      new_index = rowcol(h, cM_Size)
      if (!any(new_index %in% dist_0.1)){
        dist_0.1 = rbind(dist_0.1, new_index)
      }
    } else if(cM_dist[[j]][h] <= 0.1 & cM_dist[[j]][h] > 0.01){
      new_index = rowcol(h, cM_Size)
      if (!any(new_index %in% dist_0.01)){
        dist_0.01 = rbind(dist_0.01, new_index)
      }
    }
  }
  
  dist_100 = as.data.frame(dist_100[-1, ])
  dist_10 = as.data.frame(dist_10[-1, ])
  dist_1 = as.data.frame(dist_1[-1, ])
  dist_0.1 = as.data.frame(dist_0.1[-1, ])
  dist_0.01 = as.data.frame(dist_0.01[-1, ])
  
  # dist_100 = apply(dist_100, 2, FUN=function(x){
  #   return(rownames(marker_info_cM_parent2[marker_info_cM_parent2$CHROM==chrom_names[[j]], ])[x])
  # })
  
  cM_dist_markers[[j]] = list(dist_100 = dist_100, 
                              dist_10 = dist_10, 
                              dist_1 = dist_1, 
                              dist_0.1 = dist_0.1, 
                              dist_0.01 = dist_0.01)
}

saveRDS(cM_dist_markers, "get_markers_dist/cM_dist_markers.rds")






# archieved
# cM_dist <- dist(marker_info_cM_parent2$cM)
# format(object.size(cM_dist), units="auto")
# # [1] "4.5 Gb"
# length(cM_dist)
# # [1] 609,615,903
# max(cM_dist)
# # [1] 6657.822
# sum(cM_dist > 1000)
# # [1] 437,236,746
# sum(cM_dist <= 1000 & cM_dist > 100)
# # [1] 152,607,772
# sum(cM_dist <= 100 & cM_dist > 10)
# # [1] 17,557,173
# sum(cM_dist <= 10 & cM_dist > 1)
# # [1] 1,989,584
# sum(cM_dist <= 1 & cM_dist > 0.1)
# # [1] 202,479
# sum(cM_dist <= 0.1 & cM_dist > 0.01)
# # [1] 18,913

# cM_matx <- as.matrix(cM_dist)
# dim(cM_matx)
# # [1] 34918 34918
# format(object.size(cM_matx), units="auto")
# # [1] "9.1 Gb"
# rownames(cM_matx) <- rownames(marker_info_cM_parent)
# colnames(cM_matx) <- rownames(marker_info_cM_parent)
# 
# cM_Size <- attr(cM_dist, "Size")
# dist_1k <- matrix(NA, nrow=1, ncol=2)
# for (i in 1:length(cM_dist)){
#   if (cM_dist[i] > 1000){
#     new_index = rowcol(i, cM_Size)
#     if (!any(new_index %in% dist_1k)){
#       dist_1k = rbind(dist_1k, new_index)
#     }
#   }
# }
# 
# 
# 
# saveRDS(dist_1k, "codes/get_markers_dist_dist_1k.rds")
