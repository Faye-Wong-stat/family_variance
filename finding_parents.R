setwd("~/family_variance/")

pedigree_genotyped <- read.table("organized_data/pedigree_genotyped.txt")
marker <- read.table("organized_data/marker.txt")
marker_info <- read.table("organized_data/marker_info.txt")

# select for validation set
# select individuals that are not parents 
pedigree_valid <- pedigree_genotyped[!(pedigree_genotyped$ID %in% 
                                         union(pedigree_genotyped$P1, pedigree_genotyped$P2)), ]
# select individuals whose parents (both) are genotyped 
pedigree_valid <- pedigree_valid[(pedigree_valid$P1 %in% pedigree_genotyped$ID) & 
                                   (pedigree_valid$P2 %in% pedigree_genotyped$ID), ]
dim(pedigree_valid)
# [1] 593   3
pedigree_parent <- pedigree_genotyped[!(pedigree_genotyped$ID %in% pedigree_valid$ID), ]
dim(pedigree_parent)



sum(pedigree_valid$ID %in% union(pedigree_valid$P1, pedigree_valid$P2))
sum(pedigree_parent$ID %in% union(pedigree_parent$P1, pedigree_parent$P2))



parents <- pedigree_parent$ID[pedigree_parent$ID %in% union(pedigree_valid$P1, pedigree_valid$P2)]
non_parents_pedi <- pedigree_parent[!pedigree_parent$ID %in% parents, ]
non_parents_pedi_C <- pedigree_parent[pedigree_parent$ID %in% parents, ]

sum(non_parents_pedi_C$ID %in% union(non_parents_pedi_C$P1, non_parents_pedi_C$P2))
non_parents_pedi_C$ID[non_parents_pedi_C$ID %in% union(non_parents_pedi_C$P1, non_parents_pedi_C$P2)]
non_parents_pedi_C[non_parents_pedi_C$ID == "00C259P002" | 
                     non_parents_pedi_C$P1 == "00C259P002" | 
                     non_parents_pedi_C$P2 == "00C259P002", ]

