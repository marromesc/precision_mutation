somaticInteractions <- function(mutMat, returnAll=T) {
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
  
  sigPairsTbl[, q.value := p.adjust(pValue, method = 'fdr')]
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

plotSomaticInteraction <- function(geneMatrix, fishertest, remove_empty = F){
  corMat <- cor(t(geneMatrix), method = c("pearson", "kendall", "spearman")[1])*0
  
  for (i in 1:nrow(corMat)) {
    for (j in 1:ncol(corMat)) {
      Ind <- which((fishertest$gene1==rownames(corMat)[i] & fishertest$gene2==colnames(corMat)[j]) | (fishertest$gene2==rownames(corMat)[i] & fishertest$gene1==colnames(corMat)[j]))
      if (length(Ind)==1) {
        if (fishertest$Event[Ind] == 'Co_Occurence' & fishertest$q.value[Ind] <= 0.05){
          corMat[i,j] <- -log10(fishertest[Ind,'q.value'])
        } else if (fishertest$Event[Ind] == 'Mutually_Exclusive' & fishertest$q.value[Ind] <= 0.05){
          corMat[i,j] <- log10(fishertest[Ind,'q.value'])
        }
      }
    }
  }
  
  if(isTRUE(remove_empty)){
    corMat <- corMat[,colSums(corMat) != 0]
    corMat <- corMat[rowSums(corMat) != 0,]
  }
  
  od = hclust(dist(corMat))$order
  cell_fun=function(j, i, x, y, w, h, fill) {
    if(i >= j) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }
  }
  Heatmap(corMat,       
          col=colorRamp2(c(-3,0,3), c("darkolivegreen", 'ghostwhite', "chocolate1")),
          cluster_columns = F, 
          cluster_rows = F,
          row_names_side = 'left',
          column_names_gp = grid::gpar(fontsize = 10),
          row_names_gp = grid::gpar(fontsize = 10),
          show_row_names = T,
          show_column_names = T,
          cell_fun = cell_fun, 
          rect_gp = gpar(type = "none"),
          heatmap_legend_param = list(title = "-log10(FDR P-value)", at = c(-3, 0, 3), 
                                      labels = c("-3 (Mutual exclusivity)", '0', "3 (Co-ocurrence)"), direction = "vertical"))
}


CaseControlCOME <- function(fishertest, geneMatrix, pheno){
  events <- as.data.frame(fishertest[fishertest$q.value <= 0.05,])
  events$p.value_caco <- NA
  events$control_event <- NA
  events$case_event <- NA
  for(i in 1:nrow(events)){
    event <- events$Event[i]
    gene1 <- events$gene1[i]
    gene2 <- events$gene2[i]
    if (event == 'Co_Occurence'){
      ct <- table(ifelse(colSums(geneMatrix[c(gene1,gene2),]) == 2, 'event', 'no_event'), pheno)
      events$p.value_caco[i] <- fisher.test(ct)$p.value 
      events$control_event[i] <- ct['event','control']
      events$case_event[i] <- ct['event','case']
    } else if (event == 'Mutually_Exclusive') {
      ct <- table(ifelse(colSums(geneMatrix[c(gene1,gene2),]) == 1, 'event', 'no_event'), pheno)
      events$p.value_caco[i] <- fisher.test(ct)$p.value 
      events$control_event[i] <- ct['event','control']
      events$case_event[i] <- ct['event','case']
    }
  }
  events$q.value_caco <- p.adjust(events$p.value_caco, method = 'fdr')
  events <- events[order(events$q.value_caco),]
  return(events)
}





