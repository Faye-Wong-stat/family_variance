# this script is the same as create_marker_list_real_data
# except this uses all individuals in the pedigree and that have been genotyped
# create_marker_list_real_data only uses individuals as progenies in the pedigree and that have been genotyped

setwd("~/family_variance/")



phased_marker_info <- readRDS("cleaned_data/phased_marker_info.rds")
phased_marker_matrix <- readRDS("cleaned_data/phased_marker_matrix.rds")
marker_info_cM_parent <- readRDS("cleaned_data/marker_info_cM_parent.rds")
chrom_names <- readRDS("create_marker_list_real_data/chrom_names.rds")

strawberry <- read.csv("provided_data/data_strawberry.csv")
pedigree <- strawberry[, c("Progeny", "MotherID", "FatherID")]
pedigree <- pedigree[complete.cases(pedigree), ]

all_ind_names <- union(union(pedigree$Progeny, pedigree$MotherID), pedigree$FatherID)
all_ind_names <- intersect(all_ind_names, colnames(phased_marker_matrix))

pedigree <- pedigree[pedigree$Progeny %in% all_ind_names, ]
pedigree <- pedigree[pedigree$MotherID %in% all_ind_names, ]
pedigree <- pedigree[pedigree$FatherID %in% all_ind_names, ]

prog_ind_names <- pedigree$Progeny



for (j in 1:length(chrom_names)){
  index_duplicated = duplicated(marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[j], ]$cM)
  print(c(j, any(index_duplicated)))
  
  if (sum(index_duplicated)!=0){
    marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[j], ][index_duplicated, ] = NA
    phased_marker_matrix[phased_marker_info$CHROM==chrom_names[j], ][index_duplicated, ] = NA
    phased_marker_info[phased_marker_info$CHROM==chrom_names[j], ][index_duplicated, ] = NA
  }
}
phased_marker_matrix <- phased_marker_matrix[!is.na(marker_info_cM_parent$cM), ]
phased_marker_info <- phased_marker_info[!is.na(marker_info_cM_parent$cM), ]
marker_info_cM_parent <- marker_info_cM_parent[!is.na(marker_info_cM_parent$cM), ]



marker_list <- vector(mode="list", length=length(chrom_names))
for (j in 1:length(chrom_names)){
  marker_list[[j]] = array(0, dim=c(sum(phased_marker_info$CHROM==chrom_names[j]), 
                                    2, 
                                    length(all_ind_names)))
}
for (j in 1:length(chrom_names)){
  for (i in 1:length(all_ind_names)){
    indiv_chrom = phased_marker_matrix[phased_marker_info$CHROM==chrom_names[j], all_ind_names[i]]
    indiv_chrom_A = gsub("\\|[01]", "", indiv_chrom)
    indiv_chrom_B = gsub("[01]\\|", "", indiv_chrom)
    marker_list[[j]][indiv_chrom_A=="1", 1, i] = 1
    marker_list[[j]][indiv_chrom_B=="1", 2, i] = 1
  }
}



chrom_info <- vector(mode="list", length=length(chrom_names))
chrom_table <- data.frame(start=NA, end=NA, length=NA)
for (i in 1:length(chrom_names)){
  marker_info = marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[i], c("CHROM", "cM")]
  chrom_info[[i]] = marker_info$cM - min(marker_info$cM)
  chrom_table[i, ] = c(min(marker_info$cM), 
                       max(marker_info$cM), 
                       (max(marker_info$cM) - min(marker_info$cM)))
}



saveRDS(marker_list, "create_marker_list_real_data2/phased_marker_list.rds")
saveRDS(chrom_info, "create_marker_list_real_data2/chrom_info.rds")
saveRDS(chrom_table, "create_marker_list_real_data2/chrom_table.rds")
saveRDS(all_ind_names, "create_marker_list_real_data2/all_ind_names.rds")
saveRDS(prog_ind_names, "create_marker_list_real_data2/prog_ind_names.rds")
saveRDS(pedigree, "create_marker_list_real_data2/pedigree.rds")



