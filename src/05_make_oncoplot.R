# Maria Roman Escorza - 2023 01 04

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

library(collapse)
library(ComplexHeatmap)
library(data.table)

source('./lib/oncoPlotDetails.R')
source('./lib/mutCountMatrix.R')

mutation_datapath <- './results/Filtered_Mutations_Compiled.csv'
meta_datapath <- './results/SampleSheet.csv'

eventDataFrame <- read.csv(mutation_datapath)
SampleSheet <- read.csv(meta_datapath)

GenesPanel <- readRDS('./data/GenesPanel.RDS')


# Mutation frequency ------------------------------------------------------

geneMatrix <- mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, GenesPanel, rm_non_aberrant_samples = T)

mut_count <- colSums(geneMatrix)
mut_count <- mut_count[order(-mut_count)]
genes_sorted <- names(mut_count)


# OncoPlot for cases and controls -----------------------------------------

# define cases and controls
eventDataFrame_case <- eventDataFrame[eventDataFrame$case_control == 'case',]
eventDataFrame_control <- eventDataFrame[eventDataFrame$case_control == 'control',]

for (status in c('case', 'control')){
  if (status == 'case'){
    eventDataFrame <- eventDataFrame_case
  } else if (status == 'control') {
    eventDataFrame <- eventDataFrame_control
  }
  
  # extract patients
  patients <- SampleSheet$patient_id[SampleSheet$case_control == status]
  
  # prepare matrix
  geneMatrix <- matrix("",nrow=length(patients),ncol=length(GenesPanel))
  rownames(geneMatrix) <- patients
  colnames(geneMatrix) <- GenesPanel
  ER <- rep("",length(patients));  Her2 <- rep("",length(patients));  Grade <- rep("",length(patients)); RT = rep('', length(patients)); Batch = rep('', length(patients))
  
  for (i in 1:nrow(geneMatrix)) {
    Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i])
    pat_id <- rownames(geneMatrix)[i]
    ER[i] <- SampleSheet$er[which(SampleSheet$patient_id == pat_id)]
    Her2[i] <- SampleSheet$her2[which(SampleSheet$patient_id == pat_id)]
    Grade[i] <- SampleSheet$grade[which(SampleSheet$patient_id == pat_id)]
    RT[i] <- SampleSheet$radiotherapy[which(SampleSheet$patient_id == pat_id)]
    Batch[i] <-  SampleSheet$cohort[which(SampleSheet$patient_id == pat_id)]
    for (j in 1:ncol(geneMatrix)) {
      Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i] & eventDataFrame$gene.knowngene==colnames(geneMatrix)[j])
      if (length(Ind)==1) {
        geneMatrix[i,j] <- eventDataFrame$exonicfunc.knowngene[Ind]
      }
      if (length(Ind)>1) {
        if (length(unique(eventDataFrame$patient_id[Ind]))==1) {  
          geneMatrix[i,j] <- "multi_hit"
        }
        else {
          geneMatrix[i,j] <- eventDataFrame$exonicfunc.knowngene[Ind[1]]
        }
      }
    }
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
  pdf(paste0("./results/oncoPrint_all_", status,".pdf"), height = 5)
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
                  row_order=genes_sorted,
                  top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot()),
                  right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE)),
                  heatmap_legend_param = list(title = "Alterations"),
                  bottom_annotation = ha,
                  alter_fun_is_vectorized = FALSE,
                  pct_gp = gpar(fontsize = 6),
                  column_title = paste0(unique(eventDataFrame$case_control), ' (n=', length(patients), ')')
  ))
  dev.off()
}
