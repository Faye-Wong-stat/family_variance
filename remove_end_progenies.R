setwd("family_variance/")



pedigree <- read.csv("provided_data/pedigree.csv", row.names = 1)
marker <- read.csv("provided_data/MarkerMatrix.csv", row.names = 1)



# only select the individuals that's genotyped
pedigree_genotyped <- pedigree[pedigree$Accession_ID %in% row.names(marker), ]
dim(pedigree_genotyped)
# [1] 999   4



# select for validation set
# select individuals with both parents (remove NA parents)
pedigree_valid <- pedigree_genotyped[complete.cases(pedigree_genotyped), ]
dim(pedigree_valid)
# [1] 965   3
# select individuals that are not parents 
pedigree_valid <- pedigree_valid[!(pedigree_valid$Accession_ID %in% 
                                     union(pedigree_valid$Integrated_P1, pedigree_valid$Integrated_P2)), ]
dim(pedigree_valid)
# [1] 842   3
{
# # are their parents in this data set (genotyped)?
# sum((pedigree_valid$Integrated_P1 %in% pedigree_genotyped$Accession_ID) | 
#       (pedigree_valid$Integrated_P2 %in% pedigree_genotyped$Accession_ID))
}
# select individuals whose parents (both) are genotyped 
pedigree_valid <- pedigree_valid[(pedigree_valid$Integrated_P1 %in% pedigree_genotyped$Accession_ID) & 
                                   (pedigree_valid$Integrated_P2 %in% pedigree_genotyped$Accession_ID), ]
dim(pedigree_valid)
# [1] 436   4



# select for testing set 
# testing <=> parents
# select individuals that are parents 
{
# # complement of individuals that are not parents are parents
# pedigree_parent <- pedigree_genotyped[!pedigree_genotyped$Accession_ID %in% pedigree_valid$Accession_ID, ]
# # do they actually have progenies?
# nrow(pedigree_parent)
# # 157
# sum((pedigree_parent$Accession_ID %in% pedigree_genotyped$Integrated_P1) |
#     (pedigree_parent$Accession_ID %in% pedigree_genotyped$Integrated_P2))
# # 132
# pedigree_parent[!((pedigree_parent$Accession_ID %in% pedigree_genotyped$Integrated_P1) |
#                   (pedigree_parent$Accession_ID %in% pedigree_genotyped$Integrated_P2)), ]
# #           X Accession_ID  Integrated_P1 Integrated_P2
# #   157   157   03C114P003           <NA>          <NA>
# #   764   764   12C018P001     08C037P603          <NA>
# #   766   766   12C020P001     08C043P004          <NA>
# #   767   767   12C021P001     08C132P608          <NA>
# #   794   794   12C039P601     08C138P003          <NA>
# #   1199 1199   61C016P099     36C003P001          <NA>
# #   1295 1295   65C272P105     56C124P012          <NA>
# #   1398 1398   70C008P101           <NA>          <NA>
# #   1860 1860   88C047P002     77C030P605          <NA>
# #   1890 1890   89C041P603     81C043P603          <NA>
# #   1893 1893   89C208P012     83C093P006          <NA>
# #   1894 1894   89C218P031           <NA>          <NA>
# #   2004 2004   92C338P001     91C400P100          <NA>
# #   2005 2005   92C339P002     83C049P001          <NA>
# #   2042 2042   93C226P002     92C401P100          <NA>
# #   2100 2100   96C048P613     92C033P601          <NA>
# #   2314 2314   Akashi_005           <NA>          <NA>
# #   5173 5173     Hayazaki           <NA>          <NA>
# #   5442 5442         Kaho           <NA>          <NA>
# #   5687 5687        Loran           <NA>          <NA>
# #   6395 6395 ORUS_1407-76 Unknown_FC_056          <NA>
# #   6710 6710     PI236579           <NA>          <NA>
# #   7394 7394        Roman           <NA>          <NA>
# #   7763 7763       Tarpan           <NA>          <NA>
# #   8807 8807   Yamagata_2           <NA>          <NA>
# sum(pedigree_genotyped$Integrated_P1 %in% "03C114P003")
# sum(pedigree_genotyped$Integrated_P2 %in% "03C114P003")
# # "03C114P003" is not used as a parent
# # select individuals that actually have progenies 
# pedigree_parent <- pedigree_parent[(pedigree_parent$Accession_ID %in% pedigree_genotyped$Integrated_P1) | 
#                                      (pedigree_parent$Accession_ID %in% pedigree_genotyped$Integrated_P2), ]
}
pedigree_parent <- pedigree_genotyped[!(pedigree_genotyped$Accession_ID %in% pedigree_valid$Accession_ID), ]
dim(pedigree_parent)
# [1] 563   4
{
# sum(pedigree_genotyped$Accession_ID %in% union(pedigree_valid$Integrated_P1, pedigree_valid$Integrated_P2))
# # 78
# # some of these parents' progeny are not in the validation set
}  



# save the data sets
write.csv(pedigree_valid, "remove_end_progenies/pedigree_end_progenies.csv", row.names=F)
write.csv(pedigree_parent, "remove_end_progenies/pedigree_parents.csv", row.names=F)
write.csv(pedigree_genotyped, "remove_end_progenies/pedigree_genotyped.csv", row.names=F)



# look at trios
parents_number <- data.frame(table(c(pedigree_valid$Integrated_P1, pedigree_valid$Integrated_P2)))
# some individuals have only been parents once 
parents_once <- parents_number[parents_number$Freq==1, ]
pedigree_parent[parents_once$Var1[parents_once$Var1 %in% pedigree_parent$Accession_ID], ]
# these individuals all have parent(s)
# there are no separate trios from the group? 




