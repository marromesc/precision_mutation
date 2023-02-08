# Maria Roman Escorza - 2023 01 04

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

system('mkdir ./results/cn_muts')

library(readr)

source('./lib/oncoPlotDetails.R')
source('./lib/mutCountMatrix.R')
source('./lib/somaticInteractions.R')

mutation_datapath <- './results/Filtered_Mutations_Compiled.csv'
meta_datapath <- './results/SampleSheet.csv'
meta_cn_datapath <- '/home/maria/albyn/copynumber/precision_copynumber/results/preprocess_cn/SamplesInfo.csv'

ts_datapath <- './data/TSs_GenesPanel.csv'


# Load data ---------------------------------------------------------------

eventDataFrame <- read.csv(mutation_datapath)
eventDataFrame <- eventDataFrame[!is.na(eventDataFrame$case_control),]
SampleSheet <- read.csv(meta_datapath)
gisticRegs <- read.csv(gistic_regs_datapath)
SampleSheet_CN <- as.data.frame(read_tsv(meta_cn_datapath))
SampleSheet_CN <- SampleSheet_CN[SampleSheet_CN$first_subseq_event%in%c('ipsilateral IBC', 'NoData') & SampleSheet_CN$surgery_final=='BCS',]

SampleSheet<-SampleSheet[SampleSheet$patient_id%in%SampleSheet_CN$patient_id,]
rownames(SampleSheet_CN) <- SampleSheet_CN$patient_id
SampleSheet_CN<-SampleSheet_CN[SampleSheet$patient_id,]

geneMatrix_all <- mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, GenesPanel, annotation = T)


# Load genes affected by copy number --------------------------------------

ts <- read.csv(ts_datapath)
ts <- ts[ts$freq_loss_cases >= 0.25 | ts$freq_loss_controls >= 0.25,]
genes <- unique(ts$GenesPanel)

for (gene in genes){
  message(gene)
  # Mutation frequency ------------------------------------------------------
  geneMatrix_gene <- cbind(geneMatrix_all[,gene], SampleSheet_CN[,paste0(gene,'_cn')])
  colnames(geneMatrix_gene) <- c('mutations', 'copy_number')
  geneMatrix_gene[,2] <- ifelse(geneMatrix_gene[,2]=='neutral', NA, geneMatrix_gene[,2])
  geneMatrix_gene[,1] <- ifelse(geneMatrix_gene[,1]=='0', NA, geneMatrix_gene[,1])
  
  
  # OncoPlot for cases and controls -----------------------------------------
  
  for (status in c('case', 'control')){
    # extract patients
    patients <- SampleSheet$patient_id[SampleSheet$case_control == status]
    
    # prepare matrix
    geneMatrix <- geneMatrix_gene[patients,]
    ER <- rep("",length(patients));  Her2 <- rep("",length(patients));  Grade <- rep("",length(patients)); RT = rep('', length(patients)); Batch = rep('', length(patients))
    
    for (i in 1:nrow(geneMatrix)) {
      Ind <- which(eventDataFrame$Entrez_Gene_Id==rownames(geneMatrix)[i])
      pat_id <- rownames(geneMatrix)[i]
      ER[i] <- SampleSheet$er[which(SampleSheet$patient_id == pat_id)]
      Her2[i] <- SampleSheet$her2[which(SampleSheet$patient_id == pat_id)]
      Grade[i] <- SampleSheet$grade[which(SampleSheet$patient_id == pat_id)]
      RT[i] <- SampleSheet$radiotherapy[which(SampleSheet$patient_id == pat_id)]
      Batch[i] <-  SampleSheet$cohort[which(SampleSheet$patient_id == pat_id)]
    }
    
    Grade <- ifelse(Grade == 3, 'High', ifelse(Grade == 2, 'Int', ifelse(Grade == 1, 'Low', 'Unknown')))
    ER <- ifelse(ER == 1, 'Positive', ifelse(ER == 0, 'Negative', 'Unknown'))
    RT <- ifelse(RT == 1, 'RT+', ifelse(RT == 0, 'RT-', 'Unknown'))
    Her2 <- ifelse(Her2 == 1, 'Positive', ifelse(Her2 == 0, 'Negative', 'Unknown'))
    
    # oncoplot annotation
    ha = HeatmapAnnotation(df = data.frame(ER=ER,Her2=Her2,Grade=Grade, RT=RT, Batch=Batch), col = list(
      ER = c("Positive" =  DarkerTurquoise, "Negative" = "paleturquoise","Unknown" = "grey"),
      Her2 = c("Positive" =  DarkerTurquoise, "Negative" = "paleturquoise","Unknown" = "grey"),
      Grade = c("High" = DarkerTurquoise, 'Int' = '#8fc9c9', "Low" =  "paleturquoise","Unknown" = "grey"),
      RT = c("RT+" = DarkerTurquoise,"RT-" = "paleturquoise","Unknown" = "grey"),
      Batch = c("NKI" = '#1565C0',"DUK" = "#6A1B9A","SLO" = "#B71C1C", 'Melbourne' = '#2E7032')),
      annotation_height =1,
      annotation_width =1,
      gp = gpar(col = "black",lwd=0.05)
    )
    
    # oncoplot
    pdf(paste0("./results/cn_muts/oncoPrint_",gene,"_", status,"_5.pdf"), height = 5)
    print(oncoPrint(t(geneMatrix),
                    col=col,
                    alter_fun=alter_fun,
                    show_column_names=TRUE,
                    remove_empty_columns = FALSE,
                    remove_empty_rows = FALSE,
                    row_names_gp = gpar(fontsize = 5),
                    column_names_gp = gpar(fontsize = 3),
                    show_pct=TRUE,
                    column_order=NULL,
                    top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot()),
                    right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE)),
                    heatmap_legend_param = list(title = "Alterations"),
                    bottom_annotation = ha,
                    alter_fun_is_vectorized = FALSE,
                    pct_gp = gpar(fontsize = 6),
                    column_title = paste0(status, ' (n=', length(patients), ')')
    ))
    dev.off()
  }
  
  geneMatrix[,1] <- ifelse(!is.na(geneMatrix[,1]), 1, 0 )
  geneMatrix[,2] <- ifelse(geneMatrix[,2] == 'loss', 1, 0 )
  geneMatrix[,2] <- ifelse(!is.na(geneMatrix[,2]), 1, 0 )
  geneMatrix <- apply(geneMatrix, 2, as.numeric)
  fishertest <- somaticInteractions(geneMatrix)
  write.table(fishertest, paste0('./results/cn_muts/come_fisher_',gene,'.tsv'), sep = '\t',quote = F, row.names = F)
  
}

