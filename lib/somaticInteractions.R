somaticInteractions <- function(mutMat, pvalue=c(0.05, 0.01), fontSize=0.8, showSigSymbols=TRUE,
                                returnAll=TRUE, geneOrder=NULL, countsFontSize=0.8, countsFontColor="black", 
                                colPal="BrBG", showSum=TRUE, plotPadj=T, colNC=9, nShiftSymbols=5, 
                                sigSymbolsSize=2, sigSymbolsFontSize=0.9, pvSymbols=c(46,42), limitColorBreaks=TRUE) {
  #pairwise fisher test source code borrowed from: https://www.nature.com/articles/ncomms6901
  interactions <- sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), function(j) {f<- try(fisher.test(mutMat[,i], mutMat[,j]), silent=TRUE); if(class(f)=="try-error") NA else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
  oddsRatio <- oddsGenes <- sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), function(j) {f<- try(fisher.test(mutMat[,i], mutMat[,j]), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
  rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(mutMat)
  
  sigPairs = which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
  sigPairs2 = which(x = 10^-abs(interactions) >= 1, arr.ind = TRUE)
  
  if(nrow(sigPairs) < 1){
    stop("No meaningful interactions found.")
  }
  
  sigPairs = rbind(sigPairs, sigPairs2)
  sigPairsTbl = data.table::rbindlist(
    lapply(X = seq_along(1:nrow(sigPairs)), function(i) {
      x = sigPairs[i,]
      g1 = rownames(interactions[x[1], x[2], drop = FALSE])
      g2 = colnames(interactions[x[1], x[2], drop = FALSE])
      tbl = as.data.frame(table(apply(X = mutMat[,c(g1, g2), drop = FALSE], 1, paste, collapse = "")))
      combn = data.frame(t(tbl$Freq))
      colnames(combn) = tbl$Var1
      pval = 10^-abs(interactions[x[1], x[2]])
      fest = oddsRatio[x[1], x[2]]
      d = data.table::data.table(gene1 = g1,
                                 gene2 = g2,
                                 pValue = pval, oddsRatio = fest)
      d = cbind(d, combn)
      d
    }), fill = TRUE)
  
  sigPairsTbl[, pAdj := p.adjust(pValue, method = 'fdr')]
  sigPairsTbl[is.na(sigPairsTbl)] = 0
  sigPairsTbl$Event = ifelse(test = sigPairsTbl$oddsRatio > 1, yes = "Co_Occurence", no = "Mutually_Exclusive")
  sigPairsTbl$pair = apply(X = sigPairsTbl[,.(gene1, gene2)], MARGIN = 1, FUN = function(x) paste(sort(unique(x)), collapse = ", "))
  sigPairsTbl[,event_ratio := `01`+`10`]
  sigPairsTbl[,event_ratio := paste0(`11`, '/', event_ratio)]
  sigPairsTblSig = sigPairsTbl[order(as.numeric(pValue))][!duplicated(pair)]
  
  sigPairsTblSig = sigPairsTblSig[!gene1 == gene2] #Remove diagonal elements
  
  if(!returnAll){
    sigPairsTblSig = sigPairsTblSig[pValue < min(pvalue)]
  }
  
  return(sigPairsTblSig)
  
}
