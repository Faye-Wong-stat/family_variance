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



# effective marker sizes and heritability 
effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)

# get the parents for simulating phenotypes
parents_names_unique <- unique(c(parents_names[, 1], parents_names[, 2]))
length(parents_names_unique)
# [1] 628

# Z matrix
Z <- marker_matrix[match(parents_names_unique, rownames(marker_matrix)), ]
alle_freq <- colMeans((Z+1), na.rm=T) / 2
minor_alle_freq <- ifelse(alle_freq <= 0.5, alle_freq, (1-alle_freq))
remove_marker_indices <- minor_alle_freq < 0.05
Z <- Z[, !remove_marker_indices]
dim(Z)
# [1]   628 30811

# Z matrix corresponding marker info matrix
phased_marker_info_Z <- phased_marker_info[!remove_marker_indices, ]
dim(phased_marker_info_Z)
# [1] 30811     9

# make a list of marker names in different chromosomes 
marker_names <- rownames(phased_marker_info)
marker_names <- marker_names[!remove_marker_indices]
length(marker_names)
# [1] 30811
marker_names_list <- lapply(chrom_names, FUN=function(x){
  ifelse(phased_marker_info[marker_names, ]$CHROM == x, marker_names, NA)
})
marker_names_list <- lapply(marker_names_list, FUN=function(x){
  x[!is.na(x)]
})
sum(lengths(marker_names_list))
# [1] 30811



# dimension of simulation
n <- length(parents_names_unique)
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
# [1] 628  20
# [1] 628  20
# [1] 628  20
# [1] 628  20
# [1] 628  20
# [1] 628  20



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
    epsilon[[i]][[j]] = sapply(var_epsilon[[i]][, j], FUN=function(x){rnorm(n, 0, x)})
  }
}
names(epsilon) <- h2s
epsilon[[1]][[1]][1:5, ]
# [,1]         [,2]         [,3]        [,4]          [,5]       [,6]
# [1,] -0.10837722 -0.005398764  0.093399663  0.04647064 -0.0551443264  0.6676286
# [2,]  0.03177050  0.023004731 -0.001543923  0.21665712  0.0002526344  0.3841585
# [3,] -0.14456470 -0.020729628  0.042333562 -0.15605454  0.0755824234 -0.8205577
# [4,]  0.27598540 -0.008630617 -0.201773221  0.38593716  0.0729910595 -0.5107242
# [5,]  0.05700522 -0.012344340 -0.065581334  0.26260954 -0.1308009345  0.2582703
# [,7]        [,8]         [,9]       [,10]         [,11]
# [1,] -0.09343770  0.04859545 -0.023888610  0.00540211  0.0087676868
# [2,] -0.13743788 -0.13890353 -0.138546216  0.10602228 -0.0066961901
# [3,]  0.10636487  0.01069819 -0.094400563  0.03108144 -0.0003498374
# [4,]  0.09654422  0.04006622 -0.114517701 -0.04056865 -0.0179638677
# [5,]  0.18464679 -0.08261454 -0.005804852 -0.03994385 -0.0374040007
# [,12]        [,13]       [,14]        [,15]       [,16]
# [1,] -0.0006951676  0.018282106  0.08937382 -0.028081339 -0.01064260
# [2,]  0.0436145760 -0.002065629  0.16856573  0.013012492  0.56863097
# [3,] -0.0257574321 -0.158461972 -0.01396160  0.048134583 -0.01134394
# [4,]  0.0424018237  0.057925918  0.14266060  0.009493866  0.33782303
# [5,]  0.0139066046 -0.011559209 -0.03236541 -0.082448536  0.21573906
# [,17]       [,18]       [,19]        [,20]
# [1,]  0.01118896 -0.32306293  0.18463642 -0.027716795
# [2,] -0.02574795 -0.36186021  0.17485343 -0.005641612
# [3,]  0.07277548 -0.03771308 -0.08327565  0.006058973
# [4,]  0.05800018 -0.42742512 -0.14079662 -0.003856546
# [5,] -0.03470159  0.23183770 -0.18555594  0.016124469




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
# [,1]      [,2]       [,3]       [,4]      [,5]       [,6]
# 00C082P003 -2.146469 0.8936817 -1.0285383 -0.6197036 0.5549473  0.8044697
# 00C139P003 -3.025593 0.1016169  0.2688490 -0.9180455 1.4315654  2.4480752
# 00C259P002 -3.201928 0.8783509 -0.2630548  3.1386427 0.5468687 -0.8205577
# 00C295P002 -1.945749 0.8904499 -1.0533182  0.8446939 1.6269189 -2.2864029
# 01C132P003 -2.816715 0.5572284 -1.8827509  2.9360662 0.4792906  0.3951114
# [,7]       [,8]      [,9]       [,10]      [,11]    [,12]
# 00C082P003 -0.03730896 -1.2072964 1.6870375  1.49971744 -0.2960561 1.395556
# 00C139P003 -0.29323339 -1.3345869 0.7697137  0.44287190 -0.3115200 1.136321
# 00C259P002  1.63324600 -2.1412852 0.4261878 -0.02635889 -0.8618368 2.139026
# 00C295P002 -0.05925129  0.5182163 0.8475473  1.00014332 -0.5761493 0.749897
# 01C132P003 -0.61015111 -1.8592357 0.9024551  1.11939491 -0.5955895 1.489935
# [,13]      [,14]      [,15]      [,16]     [,17]       [,18]
# 00C082P003 0.001603771 -1.1027827  1.2185514 -0.7116198 -3.070532 -0.32306293
# 00C139P003 0.324201096  1.0916723  1.4918121 -0.1323462 -3.107469 -0.03960402
# 00C259P002 0.110619361  0.9091450  0.5899931 -0.7403233 -3.197738  2.13489859
# 00C295P002 0.540927764  0.4837803 -1.5369168 -2.7647719  0.397404 -2.27778060
# 01C132P003 0.471442637  0.5496215  0.8266315 -2.1971165 -2.373150  0.70734723
# [,19]       [,20]
# 00C082P003  0.18463642 -0.54289642
# 00C139P003 -0.72731916 -1.03955978
# 00C259P002 -0.08327565  0.08150565
# 00C295P002 -1.36634049 -0.96232804
# 01C132P003  1.35926805 -0.13120996



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

















