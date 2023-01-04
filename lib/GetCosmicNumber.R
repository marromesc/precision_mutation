GetCosmicNumber <- function(x) {
  x <- sapply(x,function(z){
    if (is.na(z)) {
      return(0)
    }
    else if (z==".") {
      return(0)
    }
    else {
      z <- strsplit(z,"=")[[1]][3]
      z <- strsplit(z,",")[[1]]
      sum <- 0
      for (i in 1:length(z)) {
        foo <- strsplit(z[i],"\\(")[[1]]
        sum <- sum + as.numeric(foo[1])
      }
      return(sum)
    }
  }) 
  return(as.numeric(x))
}