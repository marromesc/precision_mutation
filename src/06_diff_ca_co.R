# Maria Roman Escorza - 2023 01 05

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

system('mkdir ./results/diff_mut')

library(readr)

source('./lib/forestPlotFisher.R')
source('./lib/diffMut.R')
source('./lib/mutCountMatrix.R')
source('./lib/addCN2Muts.R')

mutation_datapath <- './results/Filtered_Mutations_Compiled.csv'
meta_datapath <- './results/SampleSheet.csv'
gistic_regs_datapath <- '/home/maria/albyn/precision-CaseControl/data/copynumber/gistic_regs.csv'
meta_cn_datapath <- '/home/maria/albyn/precision-CaseControl/Tables/SamplesInfo_CN.csv'

eventDataFrame <- read.csv(mutation_datapath)
eventDataFrame <- eventDataFrame[!is.na(eventDataFrame$case_control),]
SampleSheet <- read.csv(meta_datapath)
gisticRegs <- read.csv(gistic_regs_datapath)
SampleSheet_CN <- as.data.frame(read_tsv(meta_cn_datapath))

GenesPanel <- readRDS('./data/GenesPanel.RDS')


# Differential analysis ---------------------------------------------------

geneMatrix <- mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, GenesPanel)
fisherTable <- diffMut(geneMatrix, SampleSheet)
write.table(fisherTable, './results/diff_mut/diff_mut_overall.csv', sep = ',', row.names = F)

#forest plot
if (length(which(fisherTable$adjPval < 0.05))>0){
  pdf('./results/diff_mut/diff_mut_overall.pdf')
  forestPlotFisher(fisherTable, SampleSheet, fdr=0.05)
  dev.off()
}


# Differential analysis RT ---------------------------------------------------

for (rt in c(0, 1)){
  geneMatrix_i <- geneMatrix[SampleSheet$patient_id[SampleSheet$radiotherapy == rt],]
  fisherTable <- diffMut(geneMatrix_i, SampleSheet[SampleSheet$radiotherapy == rt,])
  write.table(fisherTable, paste0('./results/diff_mut/diff_mut_RT_', rt, '.csv'), sep = ',', row.names = F)
  
  if (length(which(fisherTable$adjPval < 0.05))>0){
    pdf( paste0('./results/diff_mut/diff_mut_RT_', rt, '.pdf'))
    forestPlotFisher(fisherTable, SampleSheet = SampleSheet[SampleSheet$radiotherapy == rt,], fdr=0.05)
    dev.off()
  }
}


# Differential analysis Her2 ---------------------------------------------------

for (her2 in c(0, 1)){
  geneMatrix_i <- geneMatrix[SampleSheet$patient_id[SampleSheet$her2 == her2],]
  fisherTable <- diffMut(geneMatrix_i, SampleSheet[SampleSheet$her2 == her2,])  
  write.table(fisherTable, paste0('./results/diff_mut/diff_mut_HER2_', her2, '.csv'), sep = ',', row.names = F)
  
  if (length(which(fisherTable$adjPval < 0.05))>0){
    pdf( paste0('./results/diff_mut/diff_mut_Her2_', her2, '.pdf'))
    forestPlotFisher(fisherTable, SampleSheet = SampleSheet[SampleSheet$her2 == her2,], fdr=0.05)
    dev.off()
  }
}

# Differential analysis ER ---------------------------------------------------

for (er in c(0, 1)){
  geneMatrix_i <- geneMatrix[SampleSheet$patient_id[SampleSheet$er == er],]
  fisherTable <- diffMut(geneMatrix_i, SampleSheet[SampleSheet$er == er,])
  write.table(fisherTable, paste0('./results/diff_mut/diff_mut_ER_', er, '.csv'), sep = ',', row.names = F)
  
  if (length(which(fisherTable$adjPval < 0.05))>0){
    pdf( paste0('./results/diff_mut/diff_mut_ER_', er, '.pdf'))
    forestPlotFisher(fisherTable, SampleSheet = SampleSheet[SampleSheet$er == er,], fdr=0.05)
    dev.off()
  }
}
