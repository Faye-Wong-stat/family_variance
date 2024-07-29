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
var_Zalpha               <- readRDS("simulate_phenotypes/var_Zalpha.rds")
var_epsilon              <- readRDS("simulate_phenotypes/var_epsilon.rds")

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
  
  mean = sum(results[results$family %in% c(parent1, parent2) & 
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
saveRDS(results_10000, "view_correlation_gametes/results_10000.rds")
results_10000 <- readRDS("view_correlation_gametes/results_10000.rds")



# 3*5*20 = 300
correlation_df           <- matrix(NA, nrow=300, ncol=7)
correlation_df           <- as.data.frame(correlation_df)
colnames(correlation_df) <- c("h2s", "effective_marker_sizes", "trait_number",
                              "cor_mean", "cor_sd", 
                              "BV_mean_cor", "BV_sd_cor")

for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    for (k in 1:20){
      print((i-1)*5*20 + (j-1)*20 + k)
      A = results_10000[results_10000$h2s==h2s[i] & 
                          results_10000$effective_marker_sizes==effective_marker_sizes[j] & 
                          results_10000$trait_number==k, ]
      
      stand = sqrt(var_Zalpha[k, j])
                   # +var_epsilon[[i]][k, j])
      
      correlation_df[(i-1)*5*20 + (j-1)*20 + k, ] = 
        c(h2s[i], effective_marker_sizes[j], k, 
          cor(A$BV_mean, A$gametes_BV_mean), 
          cor(A$BV_sd, A$gametes_BV_sd), 
          sqrt(sum((A$BV_mean/stand - A$gametes_BV_mean/stand)^2) / 500), 
          sqrt(sum((A$BV_sd/stand - A$gametes_BV_sd/stand)^2) / 500) )
    }
  }
}
correlation_df$effective_marker_sizes <- factor(correlation_df$effective_marker_sizes, 
                                                levels=as.factor(effective_marker_sizes))
correlation_df$h2 <- paste("h^2 == ", correlation_df$h2s, sep="")
saveRDS(correlation_df, "view_correlation_gametes/correlation_df.rds")
correlation_df <- readRDS("view_correlation_gametes/correlation_df.rds")



p1 <- ggplot(correlation_df, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=BV_mean_cor)) + 
  # facet_wrap(~h2, labeller = label_parsed) + 
  xlab("number of causal loci") + 
  ylab("standardized RMSD\nof mean") + 
  # ylim(0.99, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  theme_minimal_grid(font_size=10) #+ 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())
#  
# save_plot(paste("view_correlation_gametes/plots/", "fammean_BV_cor.pdf", sep=""),
#           plot_grid(p1),
#           base_width=6.5, base_height=2.15)

p2 <- ggplot(correlation_df, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=BV_sd_cor)) + 
  # facet_wrap(~h2, labeller = label_parsed) + 
  xlab("number of causal loci") + 
  ylab("standardized RMSD\nof sd") + 
  # ylim(0.7, 1) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  theme_minimal_grid(font_size=10) #+ 
  # theme(strip.text.x = element_blank())

save_plot(paste("view_correlation_gametes/plots/", "famsd_BV_cor.pdf", sep=""),
          plot_grid(p1, p2, labels="auto", ncol=2, align="vh"),
          base_width=6.5, base_height=3.3)




# pdf(paste("view_correlation_gametes/plots/", "means.pdf", sep=""))
# plot(A$BV_mean, A$gametes_BV_mean)
# abline(0, 1)
# dev.off()









