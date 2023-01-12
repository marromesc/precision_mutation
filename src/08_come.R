# Maria Roman Escorza - 2023 01 11  

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

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

geneMatrix <- t(mutCountMatrix(patients = SampleSheet$patient_id, GenesPanel, rm_non_aberrant_samples = T))

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


# Pairwise DISCOVER test --------------------------------------------------

#mutual-exclusivity analysis
result.mutex <- pairwise.discover.test(events[subset,], alternative = 'less', fdr.method = 'DBH')
print(result.mutex, fdr.threshold=0.05)
result.mutex1 <- as.data.frame(result.mutex)

#co-ocurrence
result.mutex <- pairwise.discover.test(events[subset,], alternative = 'greater', fdr.method = 'DBH')
print(result.mutex, fdr.threshold=0.05)
result.mutex2 <- as.data.frame(result.mutex)

#heatmap

corMat <- cor(t(geneMatrix[subset,]), method = c("pearson", "kendall", "spearman")[1])*0
for (i in 1:nrow(corMat)) {
  for (j in 1:ncol(corMat)) {
    Ind <- which((result.mutex1$gene1==rownames(corMat)[i] & result.mutex1$gene2==colnames(corMat)[j]) | (result.mutex1$gene2==rownames(corMat)[i] & result.mutex1$gene1==colnames(corMat)[j]))
    if (length(Ind)==1) {
      #corMat[i,j] <- -1 #mutual exclusivity
      corMat[i,j] <- log10(result.mutex1[Ind,'q.value'])
    }
    
    Ind <- which((result.mutex2$gene1==rownames(corMat)[i] & result.mutex2$gene2==colnames(corMat)[j]) | (result.mutex2$gene2==rownames(corMat)[i] & result.mutex2$gene1==colnames(corMat)[j]))
    if (length(Ind)==1) {
      #corMat[i,j] <- 1 #co-ocurrence
      corMat[i,j] <- -log10(result.mutex2[Ind,'q.value'])
    }
  }
}

#removing no co-ocurrent and no mutually exclusive
#corMat <- corMat[,colSums(corMat) != 0]
#corMat <- corMat[rowSums(corMat) != 0,]

pdf('./results/come_DISCOVER.pdf', width = 10, height = 10)
Heatmap(corMat,       
        col=colorRamp2(c(-3,0,3), c("darkolivegreen", 'white', "chocolate1")),
        cluster_columns = F, 
        cluster_rows = F,
        #width = ncol(corMat)*unit(4, "mm"), 
        #height = nrow(corMat)*unit(4, "mm"),
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 10),
        show_row_names = T,
        show_column_names = T,
        rect_gp = gpar(col = "grey", lwd = 1),
        heatmap_legend_param = list(title = "", at = c(-3, 3), 
                                    labels = c("Mutual exclusivity", "Co-ocurrence"), direction = "vertical")
)
dev.off()


# Pairwise fisher test ----------------------------------------------------

mutMat <- t(geneMatrix[subset,])

pdf('./results/come_fisher.pdf')
fishertest <- somaticInteractions(mutMat)
dev.off()

write.table(fishertest, './results/come_fisher_test.tsv', sep = '\t',quote = F, row.names = F)



