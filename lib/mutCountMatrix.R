mutCountMatrix <- function(patients, GenesPanel, rm_non_aberrant_samples = F){
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
  
  geneMatrix[geneMatrix > 0 ] <- 1
  
  if (isTRUE(rm_non_aberrant_samples)){
    geneMatrix <- geneMatrix[rowSums(geneMatrix) > 1,]
  }
  return(geneMatrix)
}