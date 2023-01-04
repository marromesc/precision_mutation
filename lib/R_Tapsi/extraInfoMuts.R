GetCosmicNumber <- function(x)
{
  x <- sapply(x,function(z){
    z <- strsplit(z,"cosmic70=")[[1]][2]
    z <- strsplit(z,";")[[1]][1]
    
    if (is.na(z))
    {
      return(0)
    }
    else if (z==".")
    {
      return(0)
    }
    else
    {
      z <- strsplit(z,"x3d")[[1]][3]
      z <- strsplit(z,",")[[1]]
      sum <- 0
      for (i in 1:length(z))
      {
        foo <- strsplit(z[i],"\\(")[[1]]
        sum <- sum + as.numeric(foo[1])
      }
      return(sum)
    }
  }) 
  return(as.numeric(x))
}

GetDP<- function(x)
{
  for (i in 1:length(x)){
    z = x[i]
    x[i] <- as.numeric(strsplit(strsplit(z,"DP=")[[1]][2],";")[[1]][1])
  }
  
  return(x)
}

GetFunc.refGene<- function(x)
{
  for (i in 1:length(x)){
    z = x[i]
    x[i] <- strsplit(strsplit(z,"Func.refGene=")[[1]][2],";")[[1]][1]
  }
  
  return(x)
}

GetExonicFunc.refGene<- function(x)
{
  for (i in 1:length(x)){
    z = x[i]
    x[i] <- strsplit(strsplit(z,"ExonicFunc.refGene=")[[1]][2],";")[[1]][1]
  }
  
  return(x)
}

GetGene.refGene<- function(x)
{
  for (i in 1:length(x)){
    z = x[i]
    x[i] <- strsplit(strsplit(z,"Gene.refGene=")[[1]][2],";")[[1]][1]
  }
  
  return(x)
}
