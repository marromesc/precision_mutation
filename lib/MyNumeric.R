MyNumeric <- function(x) {
  x <- sapply(x,function(z){
    if (z=="." | is.na(z) | z =='NA') {
      return(0)
    }
    else {
      return(as.numeric(z))
    }
  }) 
  return(as.numeric(x))
}