addCN2Muts <- function(gisticRegs, geneMatrix, SampleSheet_CN, onco_cn = c('amp', 'gain'), annotation=F,
                       gene_annotation = readRDS("/home/argy/Documents/rubic_run/all_genes_grch38.rds"),
                       ts=read.csv('/mnt/albyn/common/master/Human_TSGs.csv'),
                       oncogene=read.csv('/mnt/albyn/common/master/ongene_human.csv')){
  gene_annotation$Chromosome <- gsub("X", "23", gene_annotation$Chromosome)
  
  gisticRegs$gene <- "NA"
  
  genes_cn <- character()
  for (i in 1:nrow(gisticRegs)) {
    sub_we_chr <- subset(gene_annotation, gisticRegs$chr[i] == gene_annotation$Chromosome)
    subb <- subset(sub_we_chr, (sub_we_chr$Start >= gisticRegs$start[i] & sub_we_chr$Start <= gisticRegs$end[i]) | (sub_we_chr$End >= gisticRegs$start[i] & sub_we_chr$End <= gisticRegs$end[i]))
    list <- paste(c(unique(subb$Name)), collapse=", ")
    if (nrow(subb) == 0) {
      list <- "no genes"
    }
    gisticRegs$gene[i] <- list
    
    genes_cn <- c(genes_cn, subb$Name[subb$Name %in% GenesPanel])
  }
  
  genes_cn <- unique(genes_cn)
  
  SampleSheet_CN <- SampleSheet_CN[SampleSheet_CN$patient_id %in% rownames(geneMatrix),]
  geneMatrixCN <- geneMatrix[SampleSheet_CN$patient_id,]
  for (gene in genes_cn){
    message(gene)
    cn_i <- SampleSheet_CN[,paste0(gene, '_cn')] 
    
    if(gene %in% ts$GeneSymbol){
      if (isTRUE(annotation)){
        cn_i <- ifelse(cn_i %in% c('gain', 'amp'), 'neutral', cn_i)
        geneMatrixCN[,gene] <- cn_i
      } else {
        cn_i <- ifelse(cn_i == 'loss', 1, 0)
        geneMatrixCN[,gene] <- geneMatrixCN[,gene] + cn_i
      }

    } else if (gene %in% oncogene$OncogeneName){
      if (isTRUE(annotation)){
        cn_i <- ifelse(cn_i %in% c('loss'), 'neutral', cn_i)
        geneMatrixCN[,gene] <- cn_i
      } else {
        cn_i <- ifelse(cn_i %in% onco_cn, 1, 0)
        geneMatrixCN[,gene] <- geneMatrixCN[,gene] + cn_i
      }
    }
  }
  
  if (isTRUE(annotation)){
    geneMatrixCN <- ifelse(geneMatrixCN=='neutral', geneMatrix[SampleSheet_CN$patient_id,], geneMatrixCN)
  } else {
    geneMatrixCN[geneMatrixCN > 0 ] <- 1
  }
    
  return(geneMatrixCN)
}

