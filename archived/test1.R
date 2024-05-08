# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) == 0) {
#   stop("At least one argument must be supplied (input file).\n", call. = FALSE)
# } else {
#   print(paste0("Arg input:  ", args))
# }
# 
# library(ggplot2)
# 
# offsprings_name <- args[1]
# trait <- args[2]
# 
# print(offsprings_name)
# print(trait)
# 
# print("end of script")



setwd("~/family_variance/")
alphas <- readRDS("simulate_phenotypes/alphas.rds")
trait_number <- 1
mode(alphas[[1]][[trait_number]])
class(alphas[[1]][[trait_number]])
typeof(alphas[[1]][[trait_number]])

set.seed(1)
A <- matrix(sample(c(0,1), 10*30899, replace = T), nrow = 10)
B <- A %*% alphas[[1]][[trait_number]]
B