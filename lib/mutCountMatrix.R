mutCountMatrix <- function(patients, GenesPanel){
  geneMatrix <- matrix(0,nrow=length(patients),ncol=length(GenesPanel))
  rownames(geneMatrix) <- patients
  colnames(geneMatrix) <- GenesPanel
  
  for (i in 1:nrow(geneMatrix)) {
    Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i])
    for (j in 1:ncol(geneMatrix)) {
      Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i] & eventDataFrame$gene.knowngene==colnames(geneMatrix)[j])
      if (length(Ind)>=1) {
        geneMatrix[i,j] <- 1
      }
    }
  }
  
  return(geneMatrix)
}