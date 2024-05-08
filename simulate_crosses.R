setwd("~/family_variance/")
library(hypred)



# load the marker list
marker_list <- readRDS("create_marker_list/phased_marker_list.rds")
marker_info_cM_parent <- readRDS("create_marker_list/marker_info_cM_parent.rds")
indiv_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")



parents_names_matrix <- matrix(NA, nrow=1, ncol=2)
set.seed(1)
for (i in 1:20000){
  if(i %% 500 == 0){
    print(i)
  }
  # randomly select two individuals to cross
  parents_names = sample(indiv_names, 2)
  # parents_names
  # # [1] "16C088P040" "12C041P601"
  if (nrow(parents_names_matrix) > 10000){
    break
  } else {
    if (!any( apply(parents_names_matrix, 1, setequal, parents_names) )){
      parents_names_matrix = rbind(parents_names_matrix, parents_names)
      offspring_list = vector(mode="list", length=(length(chrom_names) + 1))
      offspring_list[[1]] = parents_names
      
      for (j in 1:length(chrom_names)){
        marker_info = marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[j], c("CHROM", "cM")]
        chrom_length = (max(marker_info$cM) - min(marker_info$cM)) / 100
        marker_numb = nrow(marker_info)
        genome = hypredGenome(num.chr = 1, len.chr = chrom_length, num.snp = marker_numb)
        genome = hypredNewMap(genome, new.map = ((marker_info$cM - min(marker_info$cM)) / 100))

        parent1A = marker_list[[j]][, 1, which(indiv_names==parents_names[1])]
        parent1B = marker_list[[j]][, 2, which(indiv_names==parents_names[1])]
        parent2A = marker_list[[j]][, 1, which(indiv_names==parents_names[2])]
        parent2B = marker_list[[j]][, 2, which(indiv_names==parents_names[2])]

        gametes1 = matrix(nrow=marker_numb, ncol=200)
        gametes2 = matrix(nrow=marker_numb, ncol=200)
        offspring = matrix(nrow=marker_numb, ncol=200)
        for (k in 1:200){
          gametes1[, k] = hypredRecombine(genome, genomeA = parent1A, genomeB = parent1B, mutate = F,
                                          block = FALSE)
          gametes2[, k] = hypredRecombine(genome, genomeA = parent2A, genomeB = parent2B, mutate = F,
                                          block = FALSE)
          offspring[, k] = paste(as.character(gametes1[, k]), as.character(gametes2[, k]), sep="|")
        }

        offspring_list[[j+1]] = offspring
      }

      # saveRDS(offspring_list, "codes/simulate_crosses_offspring_list.rds")

      saveRDS(offspring_list,
              paste("simulate_crosses/offspring_list_",
                    parents_names[1],
                    "_",
                    parents_names[2],
                    ".rds",
                    sep=""))
    }
  }
  
}

print(i)
print(nrow(parents_names_matrix))

# duplicated families?
parents_names_matrix <- parents_names_matrix[-1, ]
for (i in 1:nrow(parents_names_matrix)){
  if (any( apply(parents_names_matrix[-i, ], 1, setequal, parents_names_matrix[i, ]) )){
    print("duplicated parents")
    print(i)
    break
  }
}

# # load the file names, offspring names and parents names 
# file_names <- list.files("simulate_crosses/", pattern=".rds")
# # check to see if there are same pairs of parents
# offspring_names <- gsub("\\.rds", "", file_names)
# offspring_names <- gsub("offspring_list_", "", offspring_names)
# parents_names <- matrix(NA, nrow=length(offspring_names), ncol=2)
# parents_names[, 1] <- gsub("_.+", "", offspring_names)
# parents_names[, 2] <- gsub(".+_", "", offspring_names)


# create a txt of all file names for simulated offsprings
# use the following only once
file_names <- list.files("simulate_crosses/", pattern=".rds")
file_names <- paste("simulate_crosses/", file_names, sep="")
write(file_names, "simulate_crosses/file_names.txt")





