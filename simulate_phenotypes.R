setwd("~/family_variance/")
library(rrBLUP)
library(BGLR)

# cal_minor_alle_freq <- function(x){
#   alt_alle = sum(x==0) + 2*sum(x==1)
#   freq = alt_alle / (2*length(x))
#   return(ifelse(freq<=0.5, freq, (1-freq)))
# }



# load the marker list
phased_marker_info <- readRDS("cleaned_data/phased_marker_info.rds")
phased_marker_matrix <- readRDS("cleaned_data/phased_marker_matrix.rds")
marker_list <- readRDS("cleaned_data/phased_marker_list.rds")
marker_info_cM_parent <- readRDS("cleaned_data/marker_info_cM_parent.rds")
indiv_names <- readRDS("cleaned_data/indiv_names.rds")
chrom_names <- readRDS("cleaned_data/chrom_names.rds")

marker_matrix <- matrix(0, nrow=nrow(phased_marker_matrix), ncol=ncol(phased_marker_matrix))
colnames(marker_matrix) <- colnames(phased_marker_matrix)
rownames(marker_matrix) <- rownames(phased_marker_matrix)
marker_matrix[phased_marker_matrix=="0|0"] <- -1
marker_matrix[phased_marker_matrix=="1|1"] <- 1
marker_matrix <- t(marker_matrix)

# get the parents names 
offspring_names <- list.files("simulated_data/", pattern=".rds")
offspring_names <- gsub("\\.rds", "", offspring_names)
offspring_names <- gsub("offspring_list_", "", offspring_names)
parents_names <- matrix(NA, nrow=length(offspring_names), ncol=2)
parents_names[, 1] <- gsub("_.+", "", offspring_names)
parents_names[, 2] <- gsub(".+_", "", offspring_names)



# get the parents for simulating phenotypes
parents_names_unique <- unique(c(parents_names[, 1], parents_names[, 2]))

# effective marker sizes and heritability 
effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)

# Z matrix
Z <- marker_matrix[match(parents_names_unique, rownames(marker_matrix)), ]
alle_freq <- colMeans((Z+1), na.rm=T) / 2
minor_alle_freq <- ifelse(alle_freq <= 0.5, alle_freq, (1-alle_freq))
remove_marker_indices <- minor_alle_freq < 0.05
Z <- Z[, !remove_marker_indices]
dim(Z)
# [1]   639 30944



# dimension of simulation
n <- length(parents_names_unique)
m <- ncol(Z)
# m <- 30944

set.seed(1)
effective_marker_indices <- vector(mode="list", length=length(effective_marker_sizes))
for (i in 1:length(effective_marker_sizes)){
  effective_marker_indices[[i]] = matrix(NA, nrow=20, ncol=effective_marker_sizes[i])
  effective_marker_indices[[i]] = apply(effective_marker_indices[[i]], 1, FUN=function(x){
    sample(2:(m-1), effective_marker_sizes[i])})
}
effective_marker_indices[[1]]
{
  #       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
  # [1,] 17402 13219  8463 12258 24907  9942   879  6527 19285 21035  1949 12632
  # [2,] 24389 26110  4051 25174 18184  8230 25062 12205  7076  9639 18965 23674
  # [3,]  4776 29144 13500 17686 28490 17514 24230 21456 14663 18544  1531 28713
  # [4,] 26754 10540 11572 14263 13904 25306 22307 11042 26955  3477 15226 16045
  #      [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
  # [1,]  4634 10619 30087 17454 22247 23879 28487 16850
  # [2,]  6520 26664 21785  5523 22722 25560 13929 13285
  # [3,] 20729 27060  9889 24599 21876 18176 26216 13838
  # [4,]  1223 19243 16911 22223 25705 30189 21324  7977
}

# Z matrix corresponding marker info matrix
phased_marker_info_Z <- phased_marker_info[!remove_marker_indices, ]
dim(phased_marker_info_Z)
# [1] 30944     9

{# marker_names_list <- vector(mode="list", length=length(chrom_names))
# names(marker_names_list) <- chrom_names
# for (j in 1:length(chrom_names)){
#   marker_names_list[[j]] = rownames(phased_marker_info_Z)[phased_marker_info_Z$CHROM==chrom_names[j]]
# }
# marker_names_all = rownames(phased_marker_info_Z)
# 
# # remove neighboring markers from effective markers 
# effective_marker_neighbor_indices <- vector(mode="list", length=length(effective_marker_sizes))
# for (i in 1:length(effective_marker_sizes)){
#   upper_marker_index = ifelse(effective_marker_indices[[i]] + 5 > m, 
#                               m, effective_marker_indices[[i]] + 5)
#   lower_marker_index = ifelse(effective_marker_indices[[i]] - 5 < 1, 
#                               1, effective_marker_indices[[i]] - 5)
#   
#   marker = apply(effective_marker_indices[[i]], 2, FUN=function(x){
#     colnames(Z)[x]})
#   upper_marker = apply(upper_marker_index, 2, FUN=function(x){
#     colnames(Z)[x]})
#   lower_marker = apply(lower_marker_index, 2, FUN=function(x){
#     colnames(Z)[x]})
#     
#   marker_chrom = apply(marker, 2, FUN=function(x){
#     as.character(phased_marker_info_Z[x, "CHROM"])})
#   upper_marker_chrom = apply(upper_marker, 2, FUN=function(x){
#     as.character(phased_marker_info_Z[x, "CHROM"])})
#   lower_marker_chrom = apply(lower_marker, 2, FUN=function(x){
#     as.character(phased_marker_info_Z[x, "CHROM"])})
#   
#   upper_marker[marker_chrom != upper_marker_chrom] = 
#     as.character(sapply(marker_names_list[marker_chrom[marker_chrom != upper_marker_chrom]], tail, 1))
#   lower_marker[marker_chrom != lower_marker_chrom] = 
#     as.character(sapply(marker_names_list[marker_chrom[marker_chrom != lower_marker_chrom]], head, 1))
#   
#   upper_marker_index = apply(upper_marker, 2, FUN=function(x){
#     match(x, marker_names_all)})
#   lower_marker_index = apply(lower_marker, 2, FUN=function(x){
#     match(x, marker_names_all)})
#   
#   effective_marker_neighbor_indices[[i]] = vector(mode="list", length=20)
#   for (j in 1:20){
#     effective_marker_neighbor_matrix = rbind(lower_marker_index[, j], upper_marker_index[, j])
#     effective_marker_neighbor_indices[[i]][[j]] = 
#       as.vector(unlist(apply(effective_marker_neighbor_matrix, 2, FUN=function(x){x[1]:x[2]})))
#     effective_marker_neighbor_indices[[i]][[j]] = unique(effective_marker_neighbor_indices[[i]][[j]])
#     effective_marker_neighbor_indices[[i]][[j]] = 
#       setdiff(effective_marker_neighbor_indices[[i]][[j]], effective_marker_indices[[i]][, j])
#   }
#   
#   {# effective_marker_neighbor_matrix = 
#   #   matrix(c(lower_marker_index, 
#   #            upper_marker_index), ncol=2)
#   # 
#   # # A = c()
#   # # for (j in 1:nrow(effective_marker_neighbor_matrix)){
#   # #   A = c(A, (effective_marker_neighbor_matrix[j, 2] - effective_marker_neighbor_matrix[j, 1]))
#   # # }
#   # 
#   # # print(i)
#   # for (j in 1:nrow(effective_marker_neighbor_matrix)){
#   #   # print(j)
#   #   effective_marker_neighbor_indices[[i]] = c(effective_marker_neighbor_indices[[i]], 
#   #     effective_marker_neighbor_matrix[j, 1]:effective_marker_neighbor_matrix[j, 2])
#   # }
#   # 
#   # effective_marker_neighbor_indices[[i]] = unique(effective_marker_neighbor_indices[[i]])
#   # effective_marker_neighbor_indices[[i]] = 
#   #   setdiff(effective_marker_neighbor_indices[[i]], effective_marker_indices[[i]])
#   }
# }
# effective_marker_neighbor_indices[[1]][[1]]
# {
#   #  [1] 17397 17398 17399 17400 17401 17403 17404 17405 17406 17407 24384 24385
#   # [13] 24386 24387 24388 24390 24391 24392 24393 24394  4771  4772  4773  4774
#   # [25]  4775  4777  4778  4779  4780  4781 26749 26750 26751 26752 26753 26755
#   # [37] 26756 26757 26758 26759
# }
# 
# {# for (i in 1:length(effective_marker_sizes)){
# #   print(length(effective_marker_neighbor_indices[[i]]))
# # }
# # # [1] 389
# # # [1] 1581
# # # [1] 6290
# # # [1] 24943
# # # [1] 99721
# # # [1] 400277
# }
}

# alphas 
set.seed(1)
marker_effects <- vector(mode="list", length=length(effective_marker_sizes))
for (i in 1:length(effective_marker_sizes)){
  marker_effects[[i]] = 
    matrix(rnorm(effective_marker_sizes[i]*20, 0, 1), nrow=effective_marker_sizes[i], ncol=20)
}
marker_effects[[1]]
{
  #            [,1]       [,2]       [,3]        [,4]        [,5]        [,6]
  # [1,] -0.6264538  0.3295078  0.5757814 -0.62124058 -0.01619026  0.91897737
  # [2,]  0.1836433 -0.8204684 -0.3053884 -2.21469989  0.94383621  0.78213630
  # [3,] -0.8356286  0.4874291  1.5117812  1.12493092  0.82122120  0.07456498
  # [4,]  1.5952808  0.7383247  0.3898432 -0.04493361  0.59390132 -1.98935170
  #             [,7]       [,8]        [,9]      [,10]      [,11]      [,12]
  # [1,]  0.61982575 -0.4781501  0.38767161 -0.3942900 -0.1645236 -0.6887557
  # [2,] -0.05612874  0.4179416 -0.05380504 -0.0593134 -0.2533617 -0.7074952
  # [3,] -0.15579551  1.3586796 -1.37705956  1.1000254  0.6969634  0.3645820
  # [4,] -1.47075238 -0.1027877 -0.41499456  0.7631757  0.5566632  0.7685329
  #           [,13]      [,14]      [,15]       [,16]      [,17]      [,18]
  # [1,] -0.1123462  0.3411197 -0.3672215  2.40161776 -0.7432732  0.1532533
  # [2,]  0.8811077 -1.1293631 -1.0441346 -0.03924000  0.1887923  2.1726117
  # [3,]  0.3981059  1.4330237  0.5697196  0.68973936 -1.8049586  0.4755095
  # [4,] -0.6120264  1.9803999 -0.1350546  0.02800216  1.4655549 -0.7099464
  #           [,19]        [,20]
  # [1,]  0.6107264 -0.443291873
  # [2,] -0.9340976  0.001105352
  # [3,] -1.2536334  0.074341324
  # [4,]  0.2914462 -0.589520946
}

alphas <- vector(mode="list", length=length(effective_marker_sizes))
alphas <- lapply(alphas, FUN=function(x){vector(mode="list", length=20)})
for (i in 1:length(effective_marker_sizes)){
  alphas[[i]] = lapply(alphas[[i]], FUN=function(x){rep(0, m)})
  for (j in 1:20){
    alphas[[i]][[j]][effective_marker_indices[[i]][, j]] = marker_effects[[i]][, j]
    # alphas[[i]][[j]] = alphas[[i]][[j]][-effective_marker_neighbor_indices[[i]][[j]]]
  }
  
  {
  # alphas[[i]] = vector(mode="list", length=20)
  # 
  # for (j in 1:20){
  #   alphas[[i]][[j]]
  # }
  # 
  # alpha = matrix(0, nrow=m, ncol=20)
  # alpha[effective_marker_indices[[i]], ] = marker_effects[[i]]
  # for (j in 1:20){
  #   alphas[[i]]
  # }
  # alphas[[i]] = alphas[[i]][-effective_marker_neighbor_indices[[i]]]
  }
}

{# for (i in 1:length(effective_marker_sizes)){
#   print(length(alphas[[i]]))
# }
# # [1] 30904
# # [1] 30784
# # [1] 30313
# # [1] 28498
# # [1] 22549
# # [1] 10612
}

# Zalphas
Zalphas <- vector(mode="list", length=length(effective_marker_sizes))
Zalphas <- lapply(Zalphas, FUN=function(x){vector(mode="list", length=20)})
for (i in 1:length(effective_marker_sizes)){
  for (j in 1:20){
    Zalphas[[i]][[j]] = Z %*% alphas[[i]][[j]]
    # Zalphas[[i]][[j]] = Z[, -effective_marker_neighbor_indices[[i]][[j]]] %*% alphas[[i]][[j]]
  }
}
{# for (i in 1:length(effective_marker_sizes)){
#   for (j in 1:20){
#     print(dim(Zalphas[[i]][[j]]))
#   }
# }
}
Zalphas_ls_mtx <- vector(mode="list", length=length(effective_marker_sizes))
for (i in 1:length(effective_marker_sizes)){
  Zalphas_ls_mtx[[i]] = Zalphas[[i]][[1]]
  for (j in 2:20){
    Zalphas_ls_mtx[[i]] = cbind(Zalphas_ls_mtx[[i]], Zalphas[[i]][[j]])
  }
}
{# for (i in 1:length(effective_marker_sizes)){
#   print(dim(Zalphas_ls_mtx[[i]]))
# }
}


# var(Zalpha)
var_Zalpha <- sapply(Zalphas_ls_mtx, FUN=function(x){
  apply(x, 2, var)})
{apply(Zalphas_ls_mtx[[1]], 2, var)
apply(Zalphas_ls_mtx[[2]], 2, var)
apply(Zalphas_ls_mtx[[3]], 2, var)
apply(Zalphas_ls_mtx[[4]], 2, var)
apply(Zalphas_ls_mtx[[5]], 2, var)
apply(Zalphas_ls_mtx[[6]], 2, var)
}

# var(epsilon)
var_epsilon <- vector(mode="list", length=length(h2s))
for (i in 1:length(h2s)){
  var_epsilon[[i]] = apply(var_Zalpha, 2, FUN=function(x){((1-h2s[i]) / h2s[i]) * x})
}

# epsilons 
set.seed(1)
epsilon <- vector(mode="list", length=length(h2s))
epsilon <- lapply(epsilon, FUN=function(x){vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    epsilon[[i]][[j]] = sapply(var_epsilon[[i]][, j], FUN=function(x){rnorm(n, 0, x)})
  }
}
names(epsilon) <- h2s
epsilon[[1]][[1]][1:5, ]
{
  #             [,1]        [,2]         [,3]         [,4]         [,5]        [,6]
  # [1,] -0.08958388  0.03370940  0.006120775 -0.030242113 -0.014942269  0.06688514
  # [2,]  0.02626129 -0.14663846 -0.011058424 -0.261632965  0.024517635  0.20635106
  # [3,] -0.11949621 -0.02052428 -0.157279809 -0.044154795 -0.058126725  0.11124406
  # [4,]  0.22812768  0.04298133  0.074097134  0.233158390 -0.005664505  0.29080354
  # [5,]  0.04712013 -0.03972648  0.213765410  0.004905789  0.002639277 -0.01497922
  #             [,7]        [,8]         [,9]        [,10]        [,11]       [,12]
  # [1,] -0.03527881  0.14728665  0.038565240  0.034099020  0.006430808 -0.05269650
  # [2,]  0.14073904 -0.11028792 -0.009605548  0.098648775  0.077931185 -0.05577010
  # [3,] -0.09142605 -0.01285572  0.051599332 -0.040894030 -0.021832963  0.07380532
  # [4,]  0.36754826 -0.14453247  0.021061663  0.084591733 -0.018733661  0.04601430
  # [5,] -0.04260402 -0.02914087 -0.021835243 -0.003730434  0.012326273  0.05434743
  #            [,13]       [,14]        [,15]        [,16]        [,17]       [,18]
  # [1,] -0.08018578 -0.36676822 -0.013294753 -0.077494141  0.300240965 -0.13909447
  # [2,] -0.02232925  0.48266170  0.008927174  0.190557322 -0.040250795  0.18953738
  # [3,] -0.01505207  0.35349711 -0.015313121 -0.296030323 -0.286751130 -0.18781176
  # [4,]  0.16415498  0.07283527 -0.013335061  0.007659568 -0.254499032 -0.04874363
  # [5,] -0.03787164  0.07465564  0.043133512 -0.027418025  0.003836203 -0.15460840
  #             [,19]       [,20]
  # [1,] -0.003562253 -0.06330667
  # [2,]  0.016396078  0.02571416
  # [3,]  0.003266872  0.01994117
  # [4,]  0.020410480 -0.06333458
  # [5,] -0.020204171 -0.02411683
}



# Ys
Y <- vector(mode="list", length=length(h2s))
Y <- lapply(Y, FUN=function(x){vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    Y[[i]][[j]] = Zalphas_ls_mtx[[j]] + epsilon[[i]][[j]]
    rownames(Y[[i]][[j]]) = rownames(Z)
  }
}
Y[[1]][[1]][1:5,]
{
  #                 [,1]       [,2]       [,3]       [,4]      [,5]      [,6]
  # 00C259P002 1.6893402 -0.7081469  2.3990717 -0.6964163 0.5627688 0.8490214
  # 01C134P010 2.6408140 -0.8884947  0.9249414 -3.0975734 0.7248437 1.1998934
  # 02C143P003 0.5069576 -1.4213961 -0.2417347 -2.8800953 0.6421993 2.2374368
  # 02C156P003 2.8426804 -1.0283827  0.1585520  0.6919151 0.6946616 0.6536065
  # 02C230P001 3.2881267 -1.1110905  0.5191538 -2.8759683 0.5965406 2.8933499
  #                   [,7]      [,8]       [,9]      [,10]      [,11]       [,12]
  # 00C259P002  0.02084993  1.881329 -1.6996838 -0.2754734 -0.2711543 -0.07143596
  # 01C134P010  0.35266328  2.144483 -1.7478546  0.8891017 -0.1996539 -1.53179818
  # 02C143P003 -0.88622394  2.139128 -2.0743214 -0.6873161 -0.9963814  0.05506585
  # 02C156P003  2.67005064 -1.085270  1.0642546  0.9070809  0.5137061 -0.37667612
  # 02C230P001 -0.50663426  2.122843  0.0319698 -0.3133028  0.8480676 -1.09750692
  #                 [,13]      [,14]      [,15]       [,16]      [,17]      [,18]
  # 00C259P002  0.8009219 -1.8372510 -0.5551533 -0.83447566 -0.6318245  1.6458270
  # 01C134P010 -0.3080889 -1.7386054 -0.6005424  0.16255516  0.8326423  0.5117936
  # 02C143P003 -0.2564228  0.9008733 -0.3546735 -2.07515088 -1.0300243 -1.0510115
  # 02C156P003  0.2765012 -1.3976475 -0.9900261 -0.03158043  1.0222635  0.9834590
  # 02C230P001  0.8432361 -2.1466115 -0.2962269 -0.05542018  2.5310765 -0.1546084
  #                   [,19]       [,20]
  # 00C259P002 -0.003817904  0.01214001
  # 01C134P010 -0.335064733 -0.41647236
  # 02C143P003 -0.348193939  0.46433840
  # 02C156P003  0.020154829 -1.02070072
  # 02C230P001  0.562432649  0.49462172
}



# # train the models 
# lmms <- vector(mode="list", length=length(h2s))
# names(lmms) <- h2s
# lmms <- lapply(lmms, FUN=function(x){vector(mode="list", length=length(effective_marker_sizes))})
# for (i in 1:length(h2s)){
#   names(lmms[[i]]) = effective_marker_sizes
#   lmms[[i]] = lapply(lmms[[i]], FUN=function(x){vector(mode="list", length=20)})
#   for (j in 1:length(effective_marker_sizes)){
#     for (k in 1:20){
#       lmms[[i]][[j]][[k]] = mixed.solve(Y[[i]][[j]][, k], Z[, -effective_marker_indices[[j]][, k]])
#     }
#   }
# }
# lmms2 <- vector(mode="list", length=length(h2s))
# names(lmms2) <- h2s
# lmms2 <- lapply(lmms2, FUN=function(x){vector(mode="list", length=length(effective_marker_sizes))})
# for (i in 1:length(h2s)){
#   names(lmms2[[i]]) = effective_marker_sizes
#   lmms2[[i]] = lapply(lmms2[[i]], FUN=function(x){vector(mode="list", length=20)})
#   for (j in 1:length(effective_marker_sizes)){
#     for (k in 1:20){
#       lmms2[[i]][[j]][[k]] = mixed.solve(Y[[i]][[j]][, k], Z)
#     }
#   }
# }



bayesC <- vector(mode="list", length=length(h2s))
names(bayesC) <- h2s
bayesC <- lapply(bayesC, FUN=function(x){vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  names(bayesC[[i]]) = effective_marker_sizes
  bayesC[[i]] = lapply(bayesC[[i]], FUN=function(x){vector(mode="list", length=20)})
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      bayesC[[i]][[j]][[k]] = BGLR(Y[[i]][[j]][, k], 
                                   ETA=list(list(X=Z[, -effective_marker_indices[[j]][, k]], 
                                                      model="BayesC")), 
                                   verbose=F, 
                                   saveAt="codes/simulate_phenotypes2/")
    }
  }
}



# saveRDS(lmms, "codes/simulate_phenotypes2_lmms.rds")
# # saveRDS(lmms2, "codes/simulate_phenotypes2_lmms2.rds")
# # lmms2 <- readRDS("codes/simulate_phenotypes2_lmms2.rds")
# saveRDS(lmms2, "codes/simulate_phenotypes2_lmms_full.rds")
# saveRDS(remove_marker_indices, "codes/simulate_phenotypes2_remove_marker_indices.rds")
# saveRDS(effective_marker_indices, "codes/simulate_phenotypes2_effective_marker_indices.rds")
# # saveRDS(effective_marker_neighbor_indices, 
# #         "codes/simulate_phenotypes2_effective_marker_neighbor_indices.rds")
# saveRDS(Z, "codes/simulate_phenotypes2_Z.rds")
# saveRDS(alphas, "codes/simulate_phenotypes2_alphas.rds")
# saveRDS(Zalphas, "codes/simulate_phenotypes2_Zalphas.rds")
# saveRDS(var_Zalpha, "codes/simulate_phenotypes2_var_Zalpha.rds")
# saveRDS(var_epsilon, "codes/simulate_phenotypes2_var_epsilon.rds")
# saveRDS(epsilon, "codes/simulate_phenotypes2_epsilon.rds")
# saveRDS(Y, "codes/simulate_phenotypes2_Y.rds")

saveRDS(bayesC, "codes/simulate_phenotypes2_bayesC.rds")







# # archived 
# # remove 100 neighboring markers around QTL
# # get the marker names in each chromosome
# phased_marker_info_Z <- phased_marker_info[!remove_marker_indices, ]
# 
# marker_names_list <- vector(mode="list", length=length(chrom_names))
# names(marker_names_list) <- chrom_names
# for (j in 1:length(chrom_names)){
#   marker_names_list[[j]] = rownames(phased_marker_info_Z)[phased_marker_info_Z$CHROM==chrom_names[j]]
# }
# marker_names_all = rownames(phased_marker_info_Z)
# 
# # check the correlation between these markers and their neighboring markers
# # run this after the rest of the script
# neighbors <- c(10, 50, 100, 175, 250)
# max_cor <- vector(mode="list", length=length(neighbors))
# for (i in 1:length(neighbors)){
#   marker_index_sample <- sample(2:(m-1), 1000)
#   upper_marker_index = ifelse(marker_index_sample + neighbors[i] > m, 
#                               m, marker_index_sample + neighbors[i])
#   lower_marker_index = ifelse(marker_index_sample - neighbors[i] < 1, 
#                               1, marker_index_sample - neighbors[i])
#   
#   marker = colnames(Z)[marker_index_sample]
#   upper_marker = colnames(Z)[upper_marker_index]
#   lower_marker = colnames(Z)[lower_marker_index]
#   
#   marker_chrom = as.character(phased_marker_info_Z[marker, "CHROM"])
#   upper_marker_chrom = as.character(phased_marker_info_Z[upper_marker, "CHROM"])
#   lower_marker_chrom = as.character(phased_marker_info_Z[lower_marker, "CHROM"])
#   
#   upper_marker[marker_chrom != upper_marker_chrom] = 
#     as.character(sapply(marker_names_list[marker_chrom[marker_chrom != upper_marker_chrom]], tail, 1))
#   lower_marker[marker_chrom != lower_marker_chrom] = 
#     as.character(sapply(marker_names_list[marker_chrom[marker_chrom != lower_marker_chrom]], head, 1))
#   
#   for (j in 1:1000){
#     # max_cor[[j]] = cor(Z[, marker[j]], 
#     #                    Z[, which(marker_names_all==lower_marker[j]):
#     #                        which(marker_names_all==upper_marker[j])])
#     max_cor[[i]] = c(max_cor[[i]], max(cor(Z[, marker[j]],
#                                  Z[, c(which(marker_names_all==lower_marker[j]):
#                                          (which(marker_names_all==marker[j]) - 1),
#                                        (which(marker_names_all==marker[j]) + 1):
#                                          which(marker_names_all==upper_marker[j]))]
#     )))
#   }
# }
# 
# # pdf("plots/simulate_phenotypes2/max correlation distribution.pdf")
# # hist(max_cor)
# # dev.off()
# 
# sapply(max_cor, FUN=function(x){sum(x>0.95)})










