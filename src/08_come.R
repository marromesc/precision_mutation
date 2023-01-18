# Maria Roman Escorza - 2023 01 11  

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

system('mkdir ./results/come/')

library(discover)
library(readr)
library(ComplexHeatmap)
library(circlize)

source('./lib/somaticInteractions.R')
source('./lib/mutCountMatrix.R')

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


# Set up mutation data ----------------------------------------------------

geneMatrix <- t(mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, GenesPanel))
events <- discover.matrix(geneMatrix)
subset <- rowSums(geneMatrix) > 3 # remove mutations affecting less than 3 samples
geneMatrix <- geneMatrix[subset,]

# Pairwise DISCOVER test --------------------------------------------------

#mutual-exclusivity analysis
result.mutex <- pairwise.discover.test(events[subset,], alternative = 'less', fdr.method = 'DBH')
print(result.mutex, fdr.threshold=0.05)
result.mutex1 <- as.data.frame(result.mutex)
result.mutex1$Event <- rep('Mutually_Exclusive', nrow(result.mutex1))

#co-ocurrence
result.mutex <- pairwise.discover.test(events[subset,], alternative = 'greater', fdr.method = 'DBH')
print(result.mutex, fdr.threshold=0.05)
result.mutex2 <- as.data.frame(result.mutex)
result.mutex2$Event <- rep('Co_Occurence', nrow(result.mutex2))

#merge
fishertest <- rbind(result.mutex1, result.mutex2)

#heatmap
pdf('./results/come/come_DISCOVER.pdf', width = 10, height = 10)
plotSomaticInteraction(geneMatrix, fishertest)
dev.off()
write.table(fishertest, './results/come/come_DISCOVER.tsv', sep = '\t',quote = F, row.names = F)


# Pairwise fisher test ----------------------------------------------------

fishertest <- somaticInteractions(t(geneMatrix))
pdf('./results/come/come_fisher.pdf', width = 10, height = 10)
plotSomaticInteraction(geneMatrix, as.data.frame(fishertest))
dev.off()
write.table(fishertest, './results/come/come_fisher.tsv', sep = '\t',quote = F, row.names = F)

