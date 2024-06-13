# this script contains useful functions used in this project

# check if there are any special characters in individual names 
any_spec_char <- function(spec_char, x){
  return(any(grepl(spec_char, x)))
}
# replace "_" and "-" with "" in vectors
replace_spec_char <- function(x){
  return(gsub("-", "", gsub("_", "", x)))
}

# give the row and column, extract the index of a dist object
distindex <- function(i, j, n){
  # i < j <= n
  # i:row, j:col, n:number of observations 
  return(n*(i-1) - i*(i-1)/2 + j - i)
}
# extract the row and column index of a dist object
rowcol <- function(ix, n){
  # ix: the index of an element within the dist object 
  # n: the length of the object which was the input of dist
  nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc=n-(2*n-nr+1)*nr/2+ix+nr
  return(c(nr,nc))
}

# generate locations where phasing error/recombination can happen
introduce_error <- function(error, length){
  error_number = rpois(1, error*length)
  error_location = runif(error_number, 0, length)
  error_location = error_location[order(error_location)]
  return(error_location)
}
# generate error
get_error_location <- function(error_location, chrom_info){
  # error_location: locations of error in cM
  # chrom_info: locations of marker in cM (with length of number of markers)
  # output: vector of true and false (with length of number of markers)
  # true indicates a switch between two markers, 
  # false not
  position = sapply(error_location, FUN=function(x){
    which(x<chrom_info)[1]
  })
  output = rep(F, length=length(chrom_info))
  for (i in 1:length(position)){
    output[position[i]:length(output)] = !output[position[i]:length(output)]
  }
  return(output)
}

# extract the row and column index of an element 
# given the index of it within a square matrix
rowcol <- function(ix, n){
  # ix: the index of an element with the square matrix
  # n: number of rows and columns in the matrix
  nr = ix %% n 
  if (nr == 0){
    nr = n
    nc = ix %/% n 
  } else {
    nc = ix %/% n + 1
  }
  return(c(nr, nc))
}

# create nested list with the same architecture 
# for example, len=c(2, 3)
# length(list) -> 2
# length(list[[1]]) -> 3
create.list <- function(len){
  if(length(len) == 1){
    vector("list", len)
  } else {
    lapply(1:len[1], function(...) create.list(len[-1]))
  }
}
# credit: 
# https://stackoverflow.com/questions/17567172/nested-lists-how-to-define-the-size-before-entering-data










