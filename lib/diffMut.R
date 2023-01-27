diffMut <- function(geneMatrix, SampleSheet, m1Name = 'control', m2Name = 'case'){
  # sample size
  m1.sampleSize <- length(unique(SampleSheet$patient_id[SampleSheet$case_control == m1Name]))
  m2.sampleSize <- length(unique(SampleSheet$patient_id[SampleSheet$case_control == m2Name]))
  
  # get gene summary
  ct1 <- colSums(geneMatrix[SampleSheet$patient_id[SampleSheet$case_control == m1Name],])
  ct2 <- colSums(geneMatrix[SampleSheet$patient_id[SampleSheet$case_control == m2Name],])
  
  # remove mutations affecting less than 3 samples
  m1.genes = names(ct1[ct1 > 3])
  m2.genes = names(ct2[ct2 > 3])
  uniqueGenes = unique(c(m1.genes, m2.genes))
  
  # make frequency table
  m1.gs.comGenes = as.data.frame(ct1[names(ct1) %in% uniqueGenes])
  m2.gs.comGenes = as.data.frame(ct2[names(ct2) %in% uniqueGenes])
  
  # merge groups
  m.gs.meged = merge(m1.gs.comGenes, m2.gs.comGenes, by = 0, all = TRUE)
  
  #Set missing genes to zero
  m.gs.meged[is.na(m.gs.meged)] = 0
  m.gs.meged = as.data.frame(m.gs.meged)
  
  # fisher table
  fisherTable <- lapply(seq_len(nrow(m.gs.meged)), function(i){
    gene = m.gs.meged[i, 1]
    m1Mut = m.gs.meged[i,2]
    m2Mut = m.gs.meged[i,3]
    #print(i)
    
    ft_mat = matrix(c(m1Mut, m1.sampleSize-m1Mut, m2Mut, m2.sampleSize-m2Mut),
                    byrow = TRUE, nrow = 2)
    
    xf = fisher.test(ft_mat, conf.int = TRUE, conf.level = 0.95)
    
    pval = xf$p.value
    or = xf$estimate
    ci.up = xf$conf.int[2]
    ci.low = xf$conf.int[1]
    tdat = data.table::data.table(Aberration = gene, m1Mut , m2Mut, freq1 = m1Mut/m1.sampleSize, freq2 = m2Mut/m2.sampleSize ,pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
    tdat
  })
  
  fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, fill = TRUE)
  fisherTable = fisherTable[order(pval)]
  
  # fdr
  fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]
  colnames(fisherTable)[2:5] = c(m1Name, m2Name, paste0('freq_', m1Name), paste0('freq_', m2Name))
  
  return(as.data.frame(fisherTable))
}