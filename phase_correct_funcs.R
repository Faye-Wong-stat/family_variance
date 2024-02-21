# setwd("~/Desktop/Runcie Lab/Strawberry Outcross Project/codes/")
library(abind)

convert_geno_vec <- function(vec){
  # vec (2) one genotype (0, 1) at one marker 
  # output: 0, 1, 2
  # convert 00, 01/10, 11, into 0, 1, 2
  if (setequal(vec, c(0, 0))) {
    return(0)
  } else if (setequal(vec, c(0, 1))) {
    return(1)
  } else if (setequal(vec, c(1, 1))) {
    return(2)
  } else {
    print("unvalid genotype vector input, check genotype")
    return(NA)
  }
}
convert_geno_mtx <- function(matx){
  # matx (m(marker)*2(chromosome)) 
  # genotype (m) 0, 1, 2
  # convert 00, 01/10, 11, into 0, 1, 2 at m markers 
  genotype = apply(matx, 1, convert_geno_vec)
  return(genotype)
}

gen_off_chrom <- function(matx){
  # generate 4 possible offspring chrom based on 4 parent chrom
  # input: matrix of 4 parent chrom
  # output: array of all possible offspring chom 
  arry = array(c(matx[, c(1,3)], matx[, c(1,4)], matx[, c(2,3)], matx[, c(2,4)]), dim=c(2, 2, 4))
  return(arry)
}
gen_off_geno <- function(matx){
  # generate 4 possible offspring chrom based on 4 parent chrom
  # input: matrix of 4 parent chrom
  # output: array(2*1*k) of all possible offsprings, k=number of possible offsprings
  arry = gen_off_chrom(matx)
  arry2 = array(convert_geno_mtx(arry[, , 1]), dim=c(2, 1, 1))
  for (i in 2:4){
    arry3 = array(convert_geno_mtx(arry[, , i]), dim=c(2, 1, 1))
    logcs = apply(arry2, c(3), function(x){all(x==arry3[, , 1])})
    if (!any(logcs)){
      arry2 = abind(arry2, arry3, along=3)
    }
  }
  
  return(arry2)
}

get_index_non_two_het_par_cor_off <- function(vec){
  # vec is a vector with length 3
  # vec = c(p1, p2, off)
  # it the the genotype (0, 1, or 2) of two parents and their offspring
  if (!any(is.na(vec))){
    if (vec[1]==vec[2] & vec[1]==0){
      ifelse(vec[3]==0, T, F)
    } else if (vec[1]==vec[2] & vec[1]==2){
      ifelse(vec[3]==2, T, F)
    } else if (vec[1]==vec[2] & vec[1]==1){
      return(T)
    } else if (setequal(vec[1:2], c(0, 2))){
      ifelse(vec[3]==1, T, F)
    } else if (setequal(vec[1:2], c(0, 1))){
      ifelse(vec[3]!=2, T, F)
    } else if (setequal(vec[1:2], c(1, 2))){
      ifelse(vec[3]!=0, T, F)
    } else{
      return(NA)
    }
  } else{
    return(NA)
  }
}
get_index_cor_off <- function(vec){
  # vec is a vector with length 5
  # vec = c(p1A, p1B, p2A, p2B, off)
  # it the the genotype (0, 1) of two parents and their offspring (0, 1, 2)
  if (!any(is.na(vec))){
    vec2 = c(convert_geno_vec(vec[1:2]), convert_geno_vec(vec[3:4]), vec[5])
    return(get_index_non_two_het_par_cor_off(vec2))
  } else{
    return(NA)
  }
}

is_11 <- function(vec){
  # vec (2) one parent genotype (0, 1, 2) at two markers
  # is this parent an 1/1 genotype
  if (vec[1]==vec[2] & vec[1]==1){
    return(T)
  } else {
    return(F)
  }
}

same_parent <- function(matx1, matx){
  # matx (2*2) 1 parental genotypes (0, 1)'s at two markers 
  # matx1 (2*2) the parental genotype matx is comparing to 
  logcs = all(matx1==matx)
  if (logcs){
    return(T)
  } else {
    matx = matx[, 2:1]
    return(all(matx1==matx))
  }
}
same_par_off <- function(matx1, matx){
  # matx (2*3) 1 parental genotypes (0, 1)'s and 1 offspring genotype (0, 1, 2) at two markers 
  # matx1 (2*3) the parental and offspring genotype that matx is comparing to 
  same_par1 = same_parent(matx1[, 1:2], matx[, 1:2])
  same_off = all(matx1[, 3]==matx[, 3])
  return(all(c(same_par1, same_off)))
}
same_pars_off <- function(matx1, matx){
  # matx (2*5) 2 parental genotypes (0, 1)'s and 1 offspring genotype (0, 1, 2) at two markers 
  # matx1 (2*5) the parental and offspring genotype that matx is comparing to 
  same_par1 = same_parent(matx1[, 1:2], matx[, 1:2])
  same_par2 = same_parent(matx1[, 3:4], matx[, 3:4])
  same_off = all(matx1[, 5]==matx[, 5])
  return(all(c(same_par1, same_par2, same_off)))
}

twist_11 <- function(matx){
  # matx (2*2) 1 parental genotypes (0, 1)'s at two markers 
  matx1 = matx
  matx1[2, ] = matx[2, 2:1]
  return(matx1)
}





# par_geno <- matrix(c(0,0, 0,1, 1,0, 1,1), ncol=4)
{# all_par_geno <- c()
# for (i in 1:4){
#   for (j in 1:4){
#     for (k in 1:4){
#       for (l in 1:4){
#         all_par_geno[[(((i-1)*4 + j-1)*4 + k-1)*4 + l]] = 
#           matrix(c(par_geno[, i], par_geno[, j], par_geno[, k], par_geno[, l]), ncol=4)
#       }
#     }
#   }
# }
# names(all_par_geno) <- as.character(1:256)
}
# one_par_geno <- c()
# for (i in 1:4){
#   for (j in 1:4){
#     one_par_geno[[(i-1)*4 + j]] = array(c(par_geno[, i], par_geno[, j]), dim=c(2, 2))
#   }
# }
# names(one_par_geno) <- 1:16
# remove_parent <- c()
# for (i in 1:16){
#   repeat_parent = sapply(one_par_geno, same_parent, one_par_geno[[i]])
#   if (sum(repeat_parent)>1){
#     remove_parent = c(remove_parent, which(repeat_parent)[2])
#   }
# }
# one_par_geno <- one_par_geno[-unique(remove_parent)]
# 
# two_par_geno <- c()
# for (i in 1:10){
#   for (j in 1:10){
#     two_par_geno[[(i-1)*10+j]] = cbind(one_par_geno[[i]], one_par_geno[[j]])
#   }
# }
# names(two_par_geno) <- 1:100
# 
# par_off_arry <- c()
# for (i in 1:100){
#   four_off = gen_off_geno(two_par_geno[[i]])
#   par1_geno = convert_geno_mtx(two_par_geno[[i]][, 1:2])
#   par2_geno = convert_geno_mtx(two_par_geno[[i]][, 3:4])
#   parent_matrix = cbind(cbind(par1_geno, par2_geno), two_par_geno[[i]])
#   parent = array(parent_matrix, dim=c(2, 6, dim(four_off)[3]))
#   par_off = abind(parent, four_off, along=2)
#   par_off_arry = abind(par_off_arry, par_off, along=3)
# }






# off_geno <- c(0, 1, 2)
# all_off_geno <- c()
# for (i in 1:3){
#   for (j in 1:3){
#     all_off_geno[[(i-1)*3 + j]] = 
#       c(off_geno[i], off_geno[j])
#   }
# }







# archieved
{# convert_par_geno <- function(matx){
  #   # matx (2*4) parental genotypes (0, 1)'s at two markers 
  #   # lst (2) vectors of two parents unphased genotype (0, 1, 2)
  #   # splits two phased parental genotypes into two vectors of unphased genotype
  #   matx1 = matx[, 1:2]
  #   matx2 = matx[, 3:4]
  #   par1 = convert_geno_mtx(matx1)
  #   par2 = convert_geno_mtx(matx2)
  #   lst = list(par1=par1, par2=par2)
  #   return(lst)
  # }
  # any.parent.11 <- function(matx){
  #   # matx (2*4) 2 parental genotypes (0, 1)'s at two markers
  #   # is any parent a 1/1 genotype
  #   lst = convert_par_geno(matx)
  #   par1 = lst[[1]]
  #   par2 = lst[[2]]
  #   is.par1.11 = is.11(par1)
  #   is.par2.11 = is.11(par2)
  #   logcs = c(is.par1.11, is.par2.11)
  #   return(logcs)
  # }
  # parent.geno <- function(matx){
  #   # matx (2*2) 1 parental genotypes (0, 1)'s at two markers 
  #   # identify which parents is which genotype
  #   # matrix(c(0, 0, 1, 1), ncol=2) and matrix(c(1, 1, 0, 0), ncol=2) are considered type 1
  #   # matrix(c(0, 1, 0, 1), ncol=2) and matrix(c(1, 0, 1, 0), ncol=2) are considered type 2
  #   ifelse(matx[1, 1]==matx[2, 1], 1, 2)
  # }
  # sort.offspring <- function(vec){
  #   # vec (2) genotype (0, 1, 2) of offspring at one marker
  #   if (identical(vec, c(0, 0))){
  #     return(1)
  #   } else if (identical(vec, c(2, 2))){
  #     return(2)
  #   } else if (identical(vec, c(0, 2))){
  #     return(3)
  #   } else if (identical(vec, c(2, 0))){
  #     return(4)
  #   } else if (identical(vec, c(0, 1))){
  #     return(5)
  #   } else if (identical(vec, c(1, 0))){
  #     return(6)
  #   } else if (identical(vec, c(1, 2))){
  #     return(7)
  #   } else if (identical(vec, c(2, 1))){
  #     return(8)
  #   } else if (identical(vec, c(1, 1))){
  #     return(9)
  #   } else {
  #     print("unvalid genotype vector input, check genotype")
  #     return(NA)
  #   }
  # }
  # main <- function(num, matx){
  #   output = rep(F, 4)
  #   
  #   matx1 = matx[, 1:2]
  #   matx2 = matx[, 3:4]
  #   par1 = convert_geno_mtx(matx1)
  #   par2 = convert_geno_mtx(matx2)
  #   is.par1.11 = is.11(par1)
  #   is.par2.11 = is.11(par2)
  #   logcs = c(is.par1.11, is.par2.11)
  #   
  #   if (any(logcs)){
  #     
  #     if (num %in% 1:2){
  #       if (is.par1.11){
  #         output[1] = T
  #         output[2] = ifelse(parent.geno(par1)==1, T, F)
  #       }
  #       if (is.par2.11){
  #         output[3] = T
  #         output[4] = ifelse(parent.geno(par2)==1, T, F)
  #       }
  #     } else if (num %in% 3:4){
  #       if (is.par1.11){
  #         output[1] = T
  #         output[2] = ifelse(parent.geno(par1)==2, T, F)
  #       }
  #       if (is.par2.11){
  #         output[3] = T
  #         output[4] = ifelse(parent.geno(par2)==2, T, F)
  #       }
  #     } else if (num %in% 5:8){
  #       
  #     } else if (num %in% 9){
  #       
  #     } else{
  #       print("unvalid offspring category input, check off spring genotype")
  #     }
  #   }
  #   
  # }
}


