setwd("~/family_variance/")
library(rrBLUP)
library(BGLR)

# cal_minor_alle_freq <- function(x){
#   alt_alle = sum(x==0) + 2*sum(x==1)
#   freq = alt_alle / (2*length(x))
#   return(ifelse(freq<=0.5, freq, (1-freq)))
# }



# load the marker list
phased_marker_info <- readRDS("create_marker_list/phased_marker_info.rds")
phased_marker_matrix <- readRDS("create_marker_list/phased_marker_matrix.rds")
marker_list <- readRDS("create_marker_list/phased_marker_list.rds")
marker_info_cM_parent <- readRDS("create_marker_list/marker_info_cM_parent.rds")
indiv_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")

marker_matrix <- matrix(0, nrow=nrow(phased_marker_matrix), ncol=ncol(phased_marker_matrix))
colnames(marker_matrix) <- colnames(phased_marker_matrix)
rownames(marker_matrix) <- rownames(phased_marker_matrix)
marker_matrix[phased_marker_matrix=="0|0"] <- -1
marker_matrix[phased_marker_matrix=="1|1"] <- 1
marker_matrix <- t(marker_matrix)

# get the parents names 
offspring_names <- list.files("simulate_crosses/", pattern=".rds")
offspring_names <- gsub("\\.rds", "", offspring_names)
offspring_names <- gsub("offspring_list_", "", offspring_names)
parents_names <- matrix(NA, nrow=length(offspring_names), ncol=2)
parents_names[, 1] <- gsub("_.+", "", offspring_names)
parents_names[, 2] <- gsub(".+_", "", offspring_names)

# make sure one pair of parents only appear once
# for (i in 2:nrow(parents_names)){
#   if (any(apply(parents_names[-i, ], 1, FUN=function(x){setequal(x, parents_names[i, ])}))){
#     print(i)
#   }
# }
# nothing should be printed


# effective marker sizes and heritability 
effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)

# get the parents for simulating phenotypes
parents_names_unique <- unique(c(parents_names[, 1], parents_names[, 2]))
length(parents_names_unique)
# [1] 628

# Z matrix
Z <- marker_matrix
alle_freq <- colMeans((Z+1), na.rm=T) / 2
minor_alle_freq <- ifelse(alle_freq <= 0.5, alle_freq, (1-alle_freq))
remove_marker_indices <- minor_alle_freq < 0.05
Z <- Z[, !remove_marker_indices]
dim(Z)
# [1]   999 30859

# Z matrix corresponding marker info matrix
phased_marker_info_Z <- phased_marker_info[!remove_marker_indices, ]
dim(phased_marker_info_Z)
# [1] 30859     9

# make a list of marker names in different chromosomes 
marker_names <- rownames(phased_marker_info)
marker_names <- marker_names[!remove_marker_indices]
length(marker_names)
# [1] 30859
marker_names_list <- lapply(chrom_names, FUN=function(x){
  ifelse(phased_marker_info[marker_names, ]$CHROM == x, marker_names, NA)
})
marker_names_list <- lapply(marker_names_list, FUN=function(x){
  x[!is.na(x)]
})
sum(lengths(marker_names_list))
# [1] 30859



# dimension of simulation
n <- nrow(Z)
m <- ncol(Z)
p <- length(marker_names_list[[1]])



# sample locations of casual loci
set.seed(1)
effective_marker_indices <- vector(mode="list", length=length(effective_marker_sizes))
for (i in 1:length(effective_marker_sizes)){
  effective_marker_indices[[i]] = matrix(NA, nrow=20, ncol=effective_marker_sizes[i])
  effective_marker_indices[[i]] = apply(effective_marker_indices[[i]], 1, FUN=function(x){
    sample(2:(m-1), effective_marker_sizes[i])})
}
effective_marker_indices[[1]]
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

# alphas 
set.seed(1)
marker_effects <- vector(mode="list", length=length(effective_marker_sizes))
for (i in 1:length(effective_marker_sizes)){
  marker_effects[[i]] = 
    matrix(rnorm(effective_marker_sizes[i]*20, 0, 1), nrow=effective_marker_sizes[i], ncol=20)
}
marker_effects[[1]]
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

alphas <- vector(mode="list", length=length(effective_marker_sizes))
alphas <- lapply(alphas, FUN=function(x){vector(mode="list", length=20)})
for (i in 1:length(effective_marker_sizes)){
  alphas[[i]] = lapply(alphas[[i]], FUN=function(x){rep(0, m)})
  for (j in 1:20){
    alphas[[i]][[j]][effective_marker_indices[[i]][, j]] = marker_effects[[i]][, j]
  }
}



# Zalphas
Zalphas <- vector(mode="list", length=length(effective_marker_sizes))
Zalphas <- lapply(Zalphas, FUN=function(x){vector(mode="list", length=20)})
for (i in 1:length(effective_marker_sizes)){
  for (j in 1:20){
    Zalphas[[i]][[j]] = Z %*% alphas[[i]][[j]]
  }
}
# for (i in 1:length(effective_marker_sizes)){
#   for (j in 1:20){
#     print(dim(Zalphas[[i]][[j]]))
#   }
# }

Zalphas_ls_mtx <- vector(mode="list", length=length(effective_marker_sizes))
for (i in 1:length(effective_marker_sizes)){
  Zalphas_ls_mtx[[i]] = Zalphas[[i]][[1]]
  for (j in 2:20){
    Zalphas_ls_mtx[[i]] = cbind(Zalphas_ls_mtx[[i]], Zalphas[[i]][[j]])
  }
}
for (i in 1:length(effective_marker_sizes)){
  print(dim(Zalphas_ls_mtx[[i]]))
}
# [1] 999  20
# [1] 999  20
# [1] 999  20
# [1] 999  20
# [1] 999  20
# [1] 999  20



# var(Zalpha)
var_Zalpha <- sapply(Zalphas_ls_mtx, FUN=function(x){
  apply(x, 2, var)})
for (i in 1:length(Zalphas_ls_mtx)){
  print(apply(Zalphas_ls_mtx[[i]], 2, var))
}
var_Zalpha

# var(epsilon)
var_epsilon <- vector(mode="list", length=length(h2s))
for (i in 1:length(h2s)){
  var_epsilon[[i]] = apply(var_Zalpha, 2, FUN=function(x){((1-h2s[i]) / h2s[i]) * x})
}
dim(var_epsilon[[i]])
# [1] 20  6

# epsilons 
set.seed(1)
epsilon <- vector(mode="list", length=length(h2s))
epsilon <- lapply(epsilon, FUN=function(x){vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    epsilon[[i]][[j]] = sapply(var_epsilon[[i]][, j], FUN=function(x){rnorm(n, 0, sqrt(x))})
  }
}
names(epsilon) <- h2s
epsilon[[1]][[1]][1:5, ]
# [,1]        [,2]       [,3]       [,4]         [,5]         [,6]
# [1,] -0.22796696 -0.12673344  0.4087613  0.8486887  0.205794366  0.054603078
# [2,]  0.06682793  0.20627316 -0.1269851 -1.1269466 -0.005215684  0.329450847
# [3,] -0.30408581  0.20208701 -0.3606968 -0.5610995 -0.063034838 -0.232504758
# [4,]  0.58052375 -0.15825866 -0.7824313  0.3847896 -0.217468168 -0.001173787
# [5,]  0.11990810  0.03829921  0.6592802  0.2012718 -0.270165801  0.057637522
# [,7]         [,8]        [,9]       [,10]        [,11]      [,12]
# [1,]  0.08084794 -0.197746626  0.20298183  0.17722891  0.192837032 -0.2655559
# [2,] -0.27897570  0.327137511  0.18003552  0.07439838 -0.004470237  0.3578222
# [3,] -0.07757577  0.040933751  0.12289327  0.03946901 -0.087951288 -0.1609315
# [4,]  0.11723593  0.237093455 -0.69744841 -0.32401879 -0.005469401 -0.2371145
# [5,] -0.04061503 -0.007571189  0.01780596 -0.26799409  0.038588224  0.3761744
# [,13]      [,14]       [,15]      [,16]       [,17]      [,18]
# [1,] -0.148344611  0.9675692 -0.04346529  0.9938273 -0.36052918  0.4022409
# [2,]  0.015257216  0.8578723 -0.27518197  0.6254318  0.09649033 -0.3425222
# [3,]  0.007730006  0.8019192 -0.08234923 -1.9089134 -0.03239302 -0.5089561
# [4,] -0.174176436  0.1599815  0.07343479 -0.1480811 -0.16793711  0.8390493
# [5,]  0.014329959 -0.6305205  0.09611175 -0.7363177 -0.02504991 -0.0898353
# [,19]      [,20]
# [1,]  0.6765204 -0.1341933
# [2,] -0.3501695 -0.1218968
# [3,] -0.6434450  0.2133149
# [4,]  0.6990068  0.2116811
# [5,] -0.2799975  0.1032535



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
# [,1]       [,2]        [,3]       [,4]        [,5]
# Akashi005      2.0192991 -1.8434478 -0.71317668 -1.3660112  1.14963058
# BeaverSweet    1.6876402 -2.1694568 -1.45282814 -0.5506396  1.18213066
# Everbearing185 0.3478995 -0.5362377 -0.09030382 -3.2794897  0.28690005
# Hayazaki       0.6060552 -1.3840124 -2.02381948  1.9333153 -0.07866289
# K1             0.7718934 -1.6784152  0.26943693 -3.1383590  0.07976909
# [,6]       [,7]     [,8]       [,9]      [,10]
# Akashi005       1.830281734  0.7568024 1.221141 -1.1467548 -1.2577730
# BeaverSweet     2.393367527 -1.2856978 2.581909  1.5844180 -0.6567412
# Everbearing185 -1.302879082 -1.4921994 0.203930  1.1122812  0.4676682
# Hayazaki       -1.779119428  0.7370617 1.758769 -0.3097768 -2.5221965
# K1             -0.004638566 -0.8354129 1.932046  1.4759935 -1.7029960
# [,11]      [,12]     [,13]    [,14]      [,15]      [,16]
# Akashi005      -0.4152883 -0.8932738 0.9466836 2.985428  1.9376104 -0.7292890
# BeaverSweet     1.6670416 -0.2698957 0.9979392 2.875731 -0.1401274  3.0158117
# Everbearing185  1.2545134 -1.1532314 1.7866238 1.690415 -1.6962035 -1.2304119
# Hayazaki        0.4124159 -0.2371145 0.8085056 1.048477  1.1175694  0.5304204
# K1             -0.9616767  1.1447073 0.2282505 1.387339 -0.4736079 -2.4594340
# [,17]      [,18]      [,19]       [,20]
# Akashi005      -2.9087610  0.8777504  1.5786930 -0.64826760
# BeaverSweet    -2.6405338  0.3674243 -1.6357280 -0.49195272
# Everbearing185 -2.5806249 -2.6815678  0.2587276  0.06708587
# Hayazaki       -1.2506141  0.6046124  2.8548128 -0.82113176
# K1             -0.3644537 -2.1091936 -2.4677285 -0.26569706



# train the models
lmms <- vector(mode="list", length=length(h2s))
names(lmms) <- h2s
lmms <- lapply(lmms, FUN=function(x){vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  names(lmms[[i]]) = effective_marker_sizes
  lmms[[i]] = lapply(lmms[[i]], FUN=function(x){vector(mode="list", length=20)})
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      lmms[[i]][[j]][[k]] = mixed.solve(Y[[i]][[j]][, k], Z[, -effective_marker_indices[[j]][, k]])
    }
  }
}
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
                                   saveAt="simulate_phenotypes/")
    }
  }
}



saveRDS(lmms, "simulate_phenotypes/lmms.rds")
saveRDS(bayesC, "simulate_phenotypes/bayesC.rds")
# saveRDS(lmms2, "codes/simulate_phenotypes2_lmms2.rds")
saveRDS(remove_marker_indices, "simulate_phenotypes/remove_marker_indices.rds")
saveRDS(effective_marker_indices, "simulate_phenotypes/effective_marker_indices.rds")
saveRDS(Z, "simulate_phenotypes/Z.rds")
saveRDS(alphas, "simulate_phenotypes/alphas.rds")
saveRDS(Zalphas_ls_mtx, "simulate_phenotypes/Zalphas_ls_mtx.rds")
saveRDS(var_Zalpha, "simulate_phenotypes/var_Zalpha.rds")
saveRDS(var_epsilon, "simulate_phenotypes/var_epsilon.rds")
saveRDS(epsilon, "simulate_phenotypes/epsilon.rds")
saveRDS(Y, "simulate_phenotypes/Y.rds")

















