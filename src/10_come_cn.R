# Maria Roman Escorza - 2023 01 12 

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

library(discover)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

source('./lib/somaticInteractions.R')
source('./lib/addCN2Muts.R')
source('./lib/mutCountMatrix.R')
source('./lib/oncoPlotDetails.R')

mutation_datapath <- './results/Filtered_Mutations_Compiled.csv'
meta_datapath <- './results/SampleSheet.csv'
gistic_regs_datapath <- '/home/maria/albyn/precision-CaseControl/data/copynumber/gistic_regs.csv'
meta_cn_datapath <- '/home/maria/albyn/copynumber/precision_copynumber/results/preprocess_cn/SamplesInfo.csv'


# Load data ---------------------------------------------------------------

eventDataFrame <- read.csv(mutation_datapath)
eventDataFrame <- eventDataFrame[!is.na(eventDataFrame$case_control),]
SampleSheet <- read.csv(meta_datapath)
gisticRegs <- read.csv(gistic_regs_datapath)
SampleSheet_CN <- as.data.frame(read_tsv(meta_cn_datapath))

SampleSheet<-SampleSheet[SampleSheet$patient_id%in%SampleSheet_CN$patient_id,]
SampleSheet_CN<-SampleSheet_CN[SampleSheet_CN$patient_id%in%SampleSheet$patient_id,]

GenesPanel <- readRDS('./data/GenesPanel.RDS')


# Set up copy number data + mutations --------------------------------------

geneMatrix <- mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, GenesPanel)
geneMatrixCN <- t(addCN2Muts(gisticRegs=gisticRegs, geneMatrix=geneMatrix, SampleSheet_CN=SampleSheet_CN))
events <- discover.matrix(geneMatrixCN)
subset <- rowSums(geneMatrixCN) > 3 # remove mutations affecting less than 3 samples
geneMatrixCN <- geneMatrixCN[subset,]


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
pdf('./results/come/come_DISCOVER_cn.pdf', width = 10, height = 10)
plotSomaticInteraction(geneMatrixCN, fishertest)
dev.off()
write.table(fishertest, './results/come/come_DISCOVER_cn.tsv', sep = '\t',quote = F, row.names = F)


# Pairwise fisher test ----------------------------------------------------

fishertest <- somaticInteractions(mutMat=t(geneMatrixCN))
pdf('./results/come/come_fisher_cn.pdf', width = 10, height = 10)
plotSomaticInteraction(geneMatrix=geneMatrixCN, as.data.frame(fishertest))
dev.off()
write.table(fishertest, './results/come/come_fisher_cn.tsv', sep = '\t',quote = F, row.names = F)

#number of cooccurrences plot
number_cooccur <- rowSums(table(fishertest$gene1[fishertest$q.value <= 0.05 & fishertest$Event=='Co_Occurence'], fishertest$gene2[fishertest$q.value <= 0.05 & fishertest$Event=='Co_Occurence']))
number_cooccur <- data.frame(Muts = names(number_cooccur), number_cooccur = as.numeric(number_cooccur))

pdf('./results/come/number_cooccur_cn.pdf', height = 5)
number_cooccur %>% 
  ggplot(aes(fct_reorder(Muts,
                         number_cooccur), 
             number_cooccur))+
  geom_col() +
  labs(x="", y="Number of Co_Occurrences") + 
  coord_flip() + 
  theme(panel.background = element_blank(), 
        axis.line.x = element_line(colour = "grey50"),
        axis.text = element_text(size = 20))
dev.off()

#oncoplot of cooccurrences
geneMatrix <- mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, GenesPanel, annotation=T)
geneMatrixCN_annot <- t(addCN2Muts(gisticRegs=gisticRegs, geneMatrix=geneMatrix, SampleSheet_CN=SampleSheet_CN, annotation=T))
geneMatrixCN_annot <- ifelse(geneMatrixCN_annot=='0', NA, geneMatrixCN_annot)

pdf(paste0("./results/come/oncoPrint_cooccur.pdf"), height = 5)
print(oncoPrint(geneMatrixCN_annot[number_cooccur$Muts,], alter_fun = alter_fun))
dev.off()

#case control comes
rownames(SampleSheet) <- SampleSheet$patient_id
pheno <- SampleSheet[colnames(geneMatrixCN),'case_control']
caco <- CaseControlCOME(fishertest, geneMatrix = geneMatrixCN, pheno = pheno)
write.table(caco, './results/come/come_fisher_cn_caco.tsv', sep = '\t',quote = F, row.names = F)

pdf(paste0("./results/come/oncoPrint_cooccur_caco.pdf"), height = 2)
print(oncoPrint(geneMatrixCN_annot[c('CCND1', 'ATM'),], alter_fun = alter_fun,
                bottom_annotation = HeatmapAnnotation(df = data.frame(CaseControl=pheno),
                                                                      col=list(CaseControl=c('case'='green', 'control'='pink'))
                )))
dev.off()

