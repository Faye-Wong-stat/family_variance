# setwd("~/family_variance/")
# library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args))
}

offsprings_name <- args[1]
trait_number <- as.numeric(args[2])



# file_names <- read.delim("simulate_crosses/file_names.txt", header=F)
# i <- 499
# trait_number <- 1
# offsprings_name <- file_names[i,]

offsprings_name <- gsub(".rds", "", gsub("simulate_crosses/offspring_list_", "", offsprings_name))
offsprings_BV <- readRDS(paste("simulate_phenotypes_crosses/", 
                     offsprings_name, 
                     "_BV", 
                     ".rds", sep=""))
offsprings_predY_RR <- readRDS(paste("simulate_phenotypes_crosses/", 
                                     offsprings_name, 
                                     "_predY_RR", 
                                     ".rds", sep=""))
offsprings_predY <- readRDS(paste("simulate_phenotypes_crosses/",
                                     offsprings_name,
                                     "_predY",
                                     ".rds", sep=""))

effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)



set.seed(1)
select_number <- c(100, 50, 25, 5)
select_offsprings <- vector(mode="list", length=length(select_number))
for (i in 1:length(select_number)){
  select_offsprings[[i]] = sample(1:200, select_number[i])
}

result_df <- data.frame(family=rep(offsprings_name, 
                                   length(select_number)*length(h2s)*length(effective_marker_sizes)), 
                        number_offspring=rep(select_number, each=3*5),
                        h2s=as.character(rep(rep(h2s, each=5), 4)),
                        effective_marker_sizes=
                          as.character(rep(effective_marker_sizes, 4*3)),
                        trait_number=rep(trait_number, 60),
                        BV_mean=rep(NA, 60), 
                        BV_var=rep(NA, 60), 
                        predY_RR_mean=rep(NA, 60), 
                        predY_RR_var=rep(NA, 60), 
                        predY_mean=rep(NA, 60), 
                        predY_var=rep(NA, 60))

result_df$effective_marker_sizes <- 
  factor(result_df$effective_marker_sizes, levels=as.factor(effective_marker_sizes))
for (h in 1:length(select_number)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      result_df[(h-1)*3*5 + (i-1)*5 + j, 6:11] = c(mean(offsprings_BV[[i]][[j]][select_offsprings[[h]], ]), 
                                                   var(offsprings_BV[[i]][[j]][select_offsprings[[h]], ]), 
                                                   mean(offsprings_predY_RR[[i]][[j]][select_offsprings[[h]], ]), 
                                                   var(offsprings_predY_RR[[i]][[j]][select_offsprings[[h]], ]), 
                                                   mean(offsprings_predY[[i]][[j]][select_offsprings[[h]], ]), 
                                                   var(offsprings_predY[[i]][[j]][select_offsprings[[h]], ]))
    }
  }
}

saveRDS(result_df, paste("simulate_phenotypes_crosses2/",
                         offsprings_name,
                         "_result_df",
                         ".rds", sep=""))



