library(discover)

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

# Remove low frequent mutations -------------------------------------------

freq_tab <- table(eventDataFrame_all$gene.knowngene)
TopGenes = names(freq_tab[freq_tab > 5])


# Make count mutation table -----------------------------------------------

# import list of genes to GenesPanel 
# in our case we have 45 genes:
GenesPanel <- c("ARID1A","GATA3","PTEN","ATM","KMT2D","RB1","AKT1","CDH1","TP53","MAP2K4","NCOR1","NF1","ERBB2","BRCA1","RUNX1","CHEK2","BAP1","PIK3CA","FBXW7","MAP3K1","PIK3R1","KMT2C","NOTCH1","SF3B1","PBRM1","PDGFRA","CCND3","ESR1","ARID1B","EGFR","BRAF","FGFR1","MYC","CDKN2A","FGFR2","CCND1","ERBB3","MDM2","TBX3","BRCA2","IGF1R","CBFB","SMAD4","STK11","CCNE1") 

# extract patients
patients <- unique(eventDataFrame_all$precision_patient_id)

#make mutation matrix
eventDataFrame <- eventDataFrame_all
geneMatrix <- matrix(0,nrow=length(patients),ncol=length(TopGenes))
rownames(geneMatrix) <- patients
colnames(geneMatrix) <- TopGenes

for (i in 1:nrow(geneMatrix)) {
  Ind <- which(eventDataFrame$precision_patient_id==rownames(geneMatrix)[i])
  for (j in 1:ncol(geneMatrix)) {
    Ind <- which(eventDataFrame$precision_patient_id==rownames(geneMatrix)[i] & eventDataFrame$gene.knowngene==colnames(geneMatrix)[j])
    if (length(Ind)>=1) {
      geneMatrix[i,j] <- 1
    }
  }
}


# Estimation of backgroun matrix ------------------------------------------
events <- discover.matrix(t(geneMatrix))


# Pairwise tests ----------------------------------------------------------

#mutual-exclusivity analysis
result.mutex <- pairwise.discover.test(events[rownames(events) %in% GenesPanel,], alternative = 'less', fdr.method = 'DBH')
print(result.mutex, fdr.threshold=0.05)
result.mutex1 <- as.data.frame(result.mutex)

#co-ocurrence
result.mutex <- pairwise.discover.test(events[rownames(events) %in% GenesPanel,], alternative = 'greater', fdr.method = 'DBH')
print(result.mutex, fdr.threshold=1)
result.mutex2 <- as.data.frame(result.mutex)
