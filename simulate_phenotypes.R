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
effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)

# get the parents for simulating phenotypes
parents_names_unique <- unique(c(parents_names[, 1], parents_names[, 2]))
length(parents_names_unique)
# [1] 1007

# Z matrix
Z <- marker_matrix
alle_freq <- colMeans((Z+1), na.rm=T) / 2
minor_alle_freq <- ifelse(alle_freq <= 0.5, alle_freq, (1-alle_freq))
remove_marker_indices <- minor_alle_freq < 0.05
Z <- Z[, !remove_marker_indices]
dim(Z)
# [1]   1007 30899

# Z matrix corresponding marker info matrix
phased_marker_info_Z <- phased_marker_info[!remove_marker_indices, ]
dim(phased_marker_info_Z)
# [1] 30899     9

# make a list of marker names in different chromosomes 
marker_names <- rownames(phased_marker_info)
marker_names <- marker_names[!remove_marker_indices]
length(marker_names)
# [1] 30899
marker_names_list <- lapply(chrom_names, FUN=function(x){
  ifelse(phased_marker_info[marker_names, ]$CHROM == x, marker_names, NA)
})
marker_names_list <- lapply(marker_names_list, FUN=function(x){
  x[!is.na(x)]
})
sum(lengths(marker_names_list))
# [1] 30899



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
# [1] 1007  20
# [1] 1007  20
# [1] 1007  20
# [1] 1007  20
# [1] 1007  20
# [1] 1007  20



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
# [1] 20  5

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
# [,1]        [,2]        [,3]       [,4]       [,5]       [,6]
# [1,] -0.3786744 -0.52569112  0.49491401  1.3299600 -0.8187999 -0.5391856
# [2,]  0.1110074 -0.34272421  0.77411337  0.8309176 -0.2747216  0.4746961
# [3,] -0.5051149  0.27438601 -0.37474839 -1.5786185 -0.8109315  0.1039386
# [4,]  0.9643040 -0.14868231 -0.52464762  0.1145527  0.1428853 -0.8183973
# [5,]  0.1991785 -0.05948026  0.08594202 -0.8903264  1.1333448 -0.1584927
# [,7]       [,8]        [,9]       [,10]      [,11]       [,12]
# [1,] -0.32037460 -0.7582776  0.49591282 -0.53815833 0.33634171 -0.82312005
# [2,] -0.69903288  0.7619976 -0.10289295  0.02892091 0.01627543  0.34877304
# [3,] -0.22500400 -0.3762744 -0.13328569 -0.26237019 0.39007837  0.26464289
# [4,]  0.51338584  0.6299453 -0.11701819  0.14341785 0.00686086 -0.28965764
# [5,]  0.09097881  0.6099152 -0.07370258 -0.18292799 0.34575846  0.06254158
# [,13]      [,14]      [,15]       [,16]       [,17]       [,18]
# [1,] -0.1053101 -0.5048712 -0.3151645 -0.14592034  0.82008689  0.07519258
# [2,] -0.1953201 -0.6906128 -0.1583599 -0.87494455 -0.62064903 -0.76315973
# [3,]  0.4892238 -0.7174562  0.1647286  0.40762622  0.01793344 -0.82964291
# [4,] -0.1400750  0.3762430  0.4222016 -0.83261698 -0.51702502  0.88251850
# [5,]  0.5382067 -0.3841969  0.6215896  0.01632455  0.25594411 -0.35469869
# [,19]      [,20]
# [1,]  1.3448177 -0.1302142
# [2,]  0.5996885 -0.1771185
# [3,] -0.8027804 -0.3995433
# [4,] -0.1785716  0.2806827
# [5,] -0.1820598  0.0680666



# Ys
Y <- vector(mode="list", length=length(h2s))
Y <- lapply(Y, FUN=function(x){vector(mode="list", length=length(effective_marker_sizes))})
for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    Y[[i]][[j]] = Zalphas_ls_mtx[[j]] + epsilon[[i]][[j]]
    rownames(Y[[i]][[j]]) = rownames(Z)
  }
}
dim(Y[[1]][[1]])
# [1] 1007   20
Y[[1]][[1]][1:5,]
# [,1]       [,2]       [,3]     [,4]       [,5]       [,6]
# Akashi005      -2.1575985 -0.8551989  1.4605386 3.544660  0.1250363 -0.7505916
# BeaverSweet    -1.4587419 -0.7543757 -1.7032924 3.594241  1.1565912  0.4124200
# Everbearing185 -1.8912209 -0.7934465 -0.6801368 1.184705 -0.8109315  1.9564492
# Hayazaki        0.2046518  0.7475637 -1.1004290 2.995427 -0.3284010 -1.0298033
# K1             -1.3705708 -1.2852340 -0.3039012 2.494238  2.6710823  1.6194529
# [,7]       [,8]       [,9]      [,10]      [,11]
# Akashi005      -2.6228770  0.2250396  1.4579778 -0.4807180 -0.7785069
# BeaverSweet    -1.5499595  2.1632564 -0.5178875  1.1270732  0.2953535
# Everbearing185 -1.6957564 -1.1695782 -0.6020853 -0.9681056  1.6437049
# Hayazaki       -1.0134953  2.0312040  0.4573752  0.2008582 -0.7143259
# K1              0.5550091  2.0713825  1.3298391 -0.5197776 -0.8190300
# [,12]     [,13]      [,14]      [,15]     [,16]      [,17]
# Akashi005      -1.4508380 1.7859299 -0.7111277 -1.7939641 1.6332002  1.5633601
# BeaverSweet     0.7743928 1.8082661  1.5934477  0.6435266 1.5546754  0.3114165
# Everbearing185  0.7090021 1.6965983 -2.0530758  1.9136374 2.1475068 -2.3415061
# Hayazaki       -0.2098804 0.1012958  2.9386298 -0.2547115 0.9465036 -3.7875385
# K1              0.1235793 2.5417929  0.1977900  1.3263639 1.7002008 -2.8257771
# [,18]      [,19]      [,20]
# Akashi005       0.5507021  2.5984511 -0.2021020
# BeaverSweet    -0.9164131  0.5994328 -1.1344846
# Everbearing185  0.6330223 -0.5113342 -0.5457724
# Hayazaki        2.2212278  0.4643355  0.1344537
# K1              2.1401692 -2.0464195 -0.3008839



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
saveRDS(Zalphas, "simulate_phenotypes/Zalphas.rds")
saveRDS(Zalphas_ls_mtx, "simulate_phenotypes/Zalphas_ls_mtx.rds")
saveRDS(var_Zalpha, "simulate_phenotypes/var_Zalpha.rds")
saveRDS(var_epsilon, "simulate_phenotypes/var_epsilon.rds")
saveRDS(epsilon, "simulate_phenotypes/epsilon.rds")
saveRDS(Y, "simulate_phenotypes/Y.rds")

















