mutCountMatrix <- function(eventDataFrame, patients, GenesPanel, rm_non_aberrant_samples = F, annotation=F){
  geneMatrix <- matrix(0,nrow=length(patients),ncol=length(GenesPanel))
  rownames(geneMatrix) <- patients
  colnames(geneMatrix) <- GenesPanel
  
  if(isTRUE(annotation)){
    for (i in 1:nrow(geneMatrix)) {
      for (j in 1:ncol(geneMatrix)) {
        Ind <- which(eventDataFrame$Entrez_Gene_Id==rownames(geneMatrix)[i] & eventDataFrame$Hugo_Symbol==colnames(geneMatrix)[j])
        if (length(Ind)==1) {
          geneMatrix[i,j] <- eventDataFrame$Consequence[Ind]
        } else if (length(Ind)>1){
          geneMatrix[i,j] <- 'multi_hit'
        }
      }
    }
  } else {
    for (i in 1:nrow(geneMatrix)) {
      for (j in 1:ncol(geneMatrix)) {
        Ind <- which(eventDataFrame$Entrez_Gene_Id==rownames(geneMatrix)[i] & eventDataFrame$Hugo_Symbol==colnames(geneMatrix)[j])
        if (length(Ind)>=1) {
          geneMatrix[i,j] <- 1
        }
      }
    }
    geneMatrix[geneMatrix > 0 ] <- 1
  }

  
  if (isTRUE(rm_non_aberrant_samples)){
    geneMatrix <- geneMatrix[rowSums(geneMatrix) > 1,]
  }
  return(geneMatrix)
}