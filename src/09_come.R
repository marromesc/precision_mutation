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

oncogene_datapath <- '/home/maria/albyn/master/ongene_human.csv'
ts_datapath <- '/home/maria/albyn/master/Human_TSGs.csv'

eventDataFrame <- read.csv(mutation_datapath)
SampleSheet <- read.csv(meta_datapath)
gisticRegs <- read.csv(gistic_regs_datapath)
SampleSheet_CN <- as.data.frame(read_tsv(meta_cn_datapath))

GenesPanel <- readRDS('./data/GenesPanel.RDS')
oncogene <- read.csv(oncogene_datapath)
ts <- read.csv(ts_datapath)


# Set up mutation data ----------------------------------------------------

geneMatrix <- t(mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, GenesPanel))

events <- discover.matrix(geneMatrix)
subset <- rowSums(geneMatrix) > 3 # remove mutations affecting less than 3 samples

# We make a selection of genes that will be used in the pairwise co-occurrence and mutual exclusivity 
# analyses. Genes are selected if they are (1) located in a recurrently altered copy number segment and 
# included in a list of known cancer genes, or (2) included in a list of mutational driver genes.

# load cancer gene list from Bushman lab
cancer_gene <- read_tsv("http://www.bushmanlab.org/assets/doc/allonco_20130923.tsv")[['symbol']]

# load High-confidence mutational drivers list from Tamborero, D. et al. Comprehensive identification of mutational cancer driver genes across 12 tumor types. Sci Rep 3, 2650 (2013)
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fsrep02650/MediaObjects/41598_2013_BFsrep02650_MOESM2_ESM.zip", 'srep02650-s2.zip')
unzip('srep02650-s2.zip')
file.remove('srep02650-s2.zip')
mut_genes <- read_csv('srep02650-s3.csv')
high_conf_drivers <- mut_genes[mut_genes$`Putative Driver Category` == 'High Confidence Driver',][['Gene Symbol']]

# select genes
selected_genes <- unique(c(cancer_gene, high_conf_drivers))
subset[which(!(names(subset) %in% selected_genes))] <- F
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

