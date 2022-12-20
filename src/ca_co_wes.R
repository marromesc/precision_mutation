library(maftools)

setwd('/mnt/albyn/maria/precision_mutation')
system('mkdir ./results/WES')

# Load WES filtered mutations ---------------------------------------------

eventDataFrame_mutect <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds')
eventDataFrame_indel <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds')

eventDataFrame_mutect <- eventDataFrame_mutect[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                  'tumor_f', 'gene.knowngene', 'batch', 'precision_patient_id', 'tissue_pathology.y', 'event_type', 'first_subseq_event',
                                                  'er', 'her2', 'grade', 'radiotherapy')]
eventDataFrame_indel <- eventDataFrame_indel[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                'tumor_f', 'gene.knowngene', 'batch', 'precision_patient_id', 'tissue_pathology.y', 'event_type', 'first_subseq_event',
                                                'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all <- rbind(eventDataFrame_mutect, eventDataFrame_indel)

# rename mutations
eventDataFrame_all <- eventDataFrame_all[!eventDataFrame_all$exonicfunc.knowngene %in% c("synonymous SNV"),]
eventDataFrame_all$exonicfunc.knowngene[eventDataFrame_all$exonicfunc.knowngene %in% c(".")] <- "splicing"
eventDataFrame_all$exonicfunc.knowngene[eventDataFrame_all$exonicfunc.knowngene %in% c("nonsynonymous SNV")] <- "missense"
eventDataFrame_all$exonicfunc.knowngene[eventDataFrame_all$exonicfunc.knowngene %in% c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution")] <- "inframe_indel"
eventDataFrame_all$exonicfunc.knowngene[eventDataFrame_all$exonicfunc.knowngene %in% c("frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift"
eventDataFrame_all$exonicfunc.knowngene[eventDataFrame_all$exonicfunc.knowngene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_case <- eventDataFrame_all[eventDataFrame_all$case_control == 'case',]
eventDataFrame_control <- eventDataFrame_all[eventDataFrame_all$case_control == 'control',]


# Differential analysis ---------------------------------------------------

# sample size
m2.sampleSize <- length(unique(eventDataFrame_control$precision_patient_id))
m1.sampleSize <- length(unique(eventDataFrame_case$precision_patient_id))

# get gene summary
ct_case <- table(eventDataFrame_case$gene.knowngene)
ct_control <- table(eventDataFrame_control$gene.knowngene)

# remove mutations affecting less than 5 samples
m1.genes = names(ct_case[ct_case > 3])
m2.genes = names(ct_control[ct_control > 3])
uniqueGenes = unique(c(m1.genes, m2.genes))

# make frequency table
m1.gs.comGenes = as.data.frame(ct_case[names(ct_case) %in% uniqueGenes])
m2.gs.comGenes = as.data.frame(ct_control[names(ct_control) %in% uniqueGenes])

# merge groups
m.gs.meged = merge(m1.gs.comGenes, m2.gs.comGenes,
                   by = 'Var1', all = TRUE)

#Set missing genes to zero
m.gs.meged[is.na(m.gs.meged)] = 0
m.gs.meged = as.data.frame(m.gs.meged)

# fisher table
fisherTable = lapply(seq_len(nrow(m.gs.meged)), function(i){
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
  tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut , m2Mut,
                                pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
  tdat
})

fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, fill = TRUE)
fisherTable = fisherTable[order(pval)]

# fdr
fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]
colnames(fisherTable)[2:3] = c('case', 'control')

write.table(as.data.frame(fisherTable), './results/WES/diff_mut_all.csv', sep = ',')


# Differential analysis for panel genes only---------------------------------------------------

# import list of genes to GenesPanel 
# in our case we have 45 genes:
GenesPanel <- c("ARID1A","GATA3","PTEN","ATM","KMT2D","RB1","AKT1","CDH1","TP53","MAP2K4","NCOR1","NF1","ERBB2","BRCA1","RUNX1","CHEK2","BAP1","PIK3CA","FBXW7","MAP3K1","PIK3R1","KMT2C","NOTCH1","SF3B1","PBRM1","PDGFRA","CCND3","ESR1","ARID1B","EGFR","BRAF","FGFR1","MYC","CDKN2A","FGFR2","CCND1","ERBB3","MDM2","TBX3","BRCA2","IGF1R","CBFB","SMAD4","STK11","CCNE1") 

# get mut count for panel genes
m.gs.meged = m.gs.meged[m.gs.meged$Var1 %in% GenesPanel,]

# fisher table
fisherTable = lapply(seq_len(nrow(m.gs.meged)), function(i){
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
  tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut , m2Mut,
                                pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
  tdat
})

fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, fill = TRUE)
fisherTable = fisherTable[order(pval)]

# fdr
fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]
colnames(fisherTable)[2:3] = c('case', 'control')

write.table(as.data.frame(fisherTable), './results/WES/diff_mut_panel_genes.csv', sep = ',')
