# view_correlation_gametes.R
setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s                    <- c(0.8, 0.5, 0.2)
sf                     <- c(0.8, 0.99)
b                      <- qnorm(sf, 0, 1)
si                     <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)
indiv_names            <- readRDS("create_marker_list/indiv_names.rds")



file_names_results <- list.files("simulate_phenotypes_gametes/", pattern="_result_df.rds")
file_names         <- gsub("_result_df.rds", "", file_names_results)

Z                        <- readRDS("simulate_phenotypes/Z.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")

results_10000            <- readRDS("view_usefulness/results.rds")



temp <- lapply(file_names_results, FUN=function(x) {
  readRDS(paste("simulate_phenotypes_gametes/", x, sep=""))
})
for (i in 1:length(temp)){
  temp[[i]]$BV_sd       = sqrt(temp[[i]]$BV_var)
  temp[[i]]$predY_RR_sd = sqrt(temp[[i]]$predY_RR_var)
  temp[[i]]$predY_sd    = sqrt(temp[[i]]$predY_var)
  temp[[i]]             = temp[[i]][, -c(6, 8, 10)]
}

# organize the data
# in "results", each row is a family
# 300*1007=302100
results = matrix(NA, nrow=302100, ncol=10)
results <- as.data.frame(results)
for (i in 1:length(temp)){
  print(i)
  results[((i-1)*300 + 1) : ((i-1)*300 + 300), ] = temp[[i]]
}
colnames(results) <- colnames(temp[[1]])
results$effective_marker_sizes <- factor(results$effective_marker_sizes, labels = effective_marker_sizes)
saveRDS(results, "view_correlation_gametes/results.rds")
results <- readRDS("view_correlation_gametes/results.rds")





results_10000$gametes_BV_mean <- NA
results_10000$gametes_BV_sd   <- NA

for (i in 1:nrow(results_10000)){
  family   = results_10000[i, "family"]
  h2       = results_10000[i, "h2s"]
  no_QTL   = results_10000[i, "effective_marker_sizes"]
  trait_no = results_10000[i, "trait_number"]
  parent1  = strsplit(family, split="_")[[1]][1]
  parent2  = strsplit(family, split="_")[[1]][2]
  
  mean = mean(results[results$family %in% c(parent1, parent2) & 
                        results$h2s                   ==h2 & 
                        results$effective_marker_sizes==no_QTL & 
                        results$trait_number          ==trait_no, 
                      "BV_mean"])
  var = sum((results[results$family %in% c(parent1, parent2) & 
                       results$h2s                   ==h2 & 
                       results$effective_marker_sizes==no_QTL & 
                       results$trait_number          ==trait_no, 
                     "BV_sd"])^2)
  sd = sqrt(var)
  
  results_10000[i, c("gametes_BV_mean", "gametes_BV_sd")] = c(mean, sd)
}

correlation_df           <- matrix(NA, nrow=nrow(results_10000), ncol=9)
correlation_df           <- as.data.frame(correlation_df)
colnames(correlation_df) <- c("family", "parent1", "parent2", 
                              "h2s", "effective_marker_sizes", "trait_number", 
                              "BV_mean_cor", "BV_sd_cor", "BV_var_cor")


family   = results_10000[i, "family"]
h2       = results_10000[i, "h2s"]
no_QTL   = results_10000[i, "effective_marker_sizes"]
trait_no = results_10000[i, "trait_number"]
parent1  = strsplit(family, split="_")[[1]][1]
parent2  = strsplit(family, split="_")[[1]][2]

correlation_df[i, c("family","h2s", "effective_marker_sizes", "trait_number")] = 
  c(family, h2, no_QTL, trait_no)
correlation_df[i, c("parent1", "parent2")] = 
  c(parent1, parent2)

mean = mean(results[results$family %in% c(parent1, parent2) & 
                      results$h2s                   ==h2 & 
                      results$effective_marker_sizes==no_QTL & 
                      results$trait_number          ==trait_no, 
                    "BV_mean"])
var = sum((results[results$family %in% c(parent1, parent2) & 
                     results$h2s                   ==h2 & 
                     results$effective_marker_sizes==no_QTL & 
                     results$trait_number          ==trait_no, 
                   "BV_sd"])^2)
sd = sqrt(var)

results_10000[i, c("gametes_BV_mean", "gametes_BV_sd")] = c(mean, sd)

correlation_df[i, c("BV_mean_cor", "BV_sd_cor", "BV_var_cor")] = 
  c(cor(results_10000[i, "BV_mean"], mean), 
    cor(results_10000[i, "BV_mean"], ), 
    cor(results_10000[i, "BV_mean"], ))



# correlation between family and gametes mean and sd
results_BV <- results[results$h2s == h2s[1], ]

# this should be arranged in the order of indiv_names
BV_mean_family <- vector(mode)





# create a results table with usefulness
results_use <- results[rep(1:nrow(results), 2), ]
results_use$si <- rep(si, each=nrow(results))
results_use <- results_use[, c(1, 11, 2:10)]
results_use$BV_use <- results_use$BV_mean + results_use$si * results_use$BV_sd
results_use$predY_RR_use <- results_use$predY_RR_mean + results_use$si * results_use$predY_RR_sd
results_use$predY_use <- results_use$predY_mean + results_use$si * results_use$predY_sd











