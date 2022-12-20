library(collapse)
library(ComplexHeatmap)
library(data.table)

setwd('/mnt/albyn/maria/precision_mutation')
system('mkdir ./results/DCIS_Precision_Panel_NKI_OncoPrint')


# Oncoplot colors ---------------------------------------------------------

col = c(missense = "violet", multi_hit = "black", nonsense = "red", UTR_variant = "yellow", splice_site = "green", inframe_indel = "orange",frameshift="blue",splicing="yellow")
alter_fun <- list(
  missense=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["missense"], col=NA))
  },
  UTR_variant=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["UTR_variant"], col=NA))
  },
  splice_site=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["splice_site"], col=NA))
  },
  inframe_indel=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["inframe_indel"], col=NA))
  },
  frameshift=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["frameshift"], col=NA))
  },
  nonsense=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["nonsense"], col=NA))
  },
  multi_hit=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["multi_hit"], col=NA))
  },
  splicing=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["splicing"], col=NA))
  }
)

PastelOrange <- rgb(red=100/100, green=70/100, blue=28/100, alpha=1)
PastelRed <- rgb(red=100/100, green=41/100, blue=38/100, alpha=1)
PastelYellow <- rgb(red=99/100, green=99/100, blue=59/100, alpha=1)

PastelGreen1 <- rgb(red=218/255, green=241/255, blue=219/255, alpha=1)
PastelGreen2 <- rgb(red=143/255, green=214/255, blue=148/255, alpha=1)
PastelGreen3 <- rgb(red=60/255, green=165/255, blue=68/255, alpha=1)
PastelGreen4 <- rgb(red=14/255, green=37/255, blue=15/255, alpha=1)

PastelViolet1 <- rgb(red=229/255, green=204/255, blue=228/255, alpha=1)
PastelViolet2 <- rgb(red=177/255, green=102/255, blue=174/255, alpha=1)

PastelBlue1 <- rgb(red=22/255, green=232/255, blue=235/255, alpha=1)
PastelBlue2 <- rgb(red=126/255, green=164/255, blue=179/255, alpha=1)

DarkerTurquoise <- rgb(0/255,131/255,133/255)
Maroon <- rgb(128/255,0,0)


# Load NKI Panel samples --------------------------------------------------

samples <- as.data.frame(fread('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples.txt'))
samples <- samples[samples$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA),]
patients <- unique(samples$patient_id); length(patients)
samples$n_mut <- NA
rownames(samples) <- samples$patient_id


# Load NKI Panel filtered mutations ---------------------------------------------

eventDataFrame_mutect <- readRDS('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.rds')

eventDataFrame_all <- eventDataFrame_mutect[,c('Consequence', 'CHROM', 'POS', 'POS', 'REF', 'ALT',
                                              'vaf_DCIS_DNA', 'Gene.refGene', 'patient_id', 'tissue_pathology', 'first_subseq_event',
                                              'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

# define cases and controls
eventDataFrame_case <- eventDataFrame_all[eventDataFrame_all$case_control == 'case',]
eventDataFrame_control <- eventDataFrame_all[eventDataFrame_all$case_control == 'control',]


# Mutation frequency of panel ------------------------------------------------------

# import list of genes to GenesPanel 
# in our case we have 45 genes:
# GenesPanel <- c("ARID1A","GATA3","PTEN","ATM","KMT2D","RB1","AKT1","CDH1","TP53","MAP2K4","NCOR1","NF1","ERBB2","BRCA1","RUNX1","CHEK2","BAP1","PIK3CA","FBXW7","MAP3K1","PIK3R1","KMT2C","NOTCH1","SF3B1","PBRM1","PDGFRA","CCND3","ESR1","ARID1B","EGFR","BRAF","FGFR1","MYC","CDKN2A","FGFR2","CCND1","ERBB3","MDM2","TBX3","BRCA2","IGF1R","CBFB","SMAD4","STK11","CCNE1") 
# saveRDS(GenesPanel, './data/GenesPanel.RDS')
GenesPanel <- readRDS('./data/GenesPanel.RDS')

#make mutation matrix
# eventDataFrame <- eventDataFrame_all
# geneMatrix <- matrix(0,nrow=length(patients),ncol=length(GenesPanel))
# rownames(geneMatrix) <- patients
# colnames(geneMatrix) <- GenesPanel
# 
# for (i in 1:nrow(geneMatrix)) {
#   Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i])
#   for (j in 1:ncol(geneMatrix)) {
#     Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i] & eventDataFrame$Gene.refGene==colnames(geneMatrix)[j])
#     if (length(Ind)>=1) {
#       geneMatrix[i,j] <- 1
#     }
#   }
# }
# 
# mut_count <- colSums(geneMatrix)
# mut_count <- mut_count[order(-mut_count)]
# genes_sorted <- names(mut_count)
# saveRDS(genes_sorted, './data/GenesPanel_sorted.RDS')
genes_sorted <- readRDS('./data/GenesPanel_sorted.RDS')

# Make oncoprint for panel genes ----------------------------------------------------------

for (status in c('case', 'control')){
  if (status == 'case'){
    eventDataFrame <- eventDataFrame_case
  } else if (status == 'control') {
    eventDataFrame <- eventDataFrame_control
  }
  
  # extract patients
  patients <- unique(samples$patient_id[samples$case_control == status])
  
  # prepare matrix
  geneMatrix <- matrix("",nrow=length(patients),ncol=length(GenesPanel))
  rownames(geneMatrix) <- patients
  colnames(geneMatrix) <- GenesPanel
  ER <- rep("",length(patients));  Her2 <- rep("",length(patients));  Grade <- rep("",length(patients)); RT = rep('', length(patients)); Batch = rep('', length(patients))
  
  for (i in 1:nrow(geneMatrix)) {
    Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i])
    pat_id <- rownames(geneMatrix)[i]
    ER[i] <- samples$er[samples$patient_id == pat_id]
    Her2[i] <- samples$her2[samples$patient_id == pat_id]
    Grade[i] <- samples$grade[samples$patient_id == pat_id]
    RT[i] <- samples$radiotherapy[samples$patient_id == pat_id]
    Batch[i] <- 'NKI'
    for (j in 1:ncol(geneMatrix)) {
      Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i] & eventDataFrame$Gene.refGene==colnames(geneMatrix)[j])
      if (length(Ind)==1) {
        geneMatrix[i,j] <- eventDataFrame$Consequence[Ind]
      }
      if (length(Ind)>1) {
        if (length(unique(eventDataFrame$patient_id[Ind]))==1) {  
          geneMatrix[i,j] <- "multi_hit"
        }
        else {
          geneMatrix[i,j] <- eventDataFrame$Consequence[Ind[1]]
        }
      }
    }
  }
  
  geneMatrix_bin <- ifelse(geneMatrix == "", 0, 1)
  samples[rownames(geneMatrix_bin), 'n_mut'] <- rowSums(geneMatrix_bin)
  
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
    #numbers=anno_text(patients,gp=gpar(fontsize=3)),
    annotation_height =1,
    annotation_width =1,
    #annotation_name_gp = gpar(fontsize=8),
    gp = gpar(col = "black",lwd=0.05)
  )

  # oncoplot
  pdf(paste0("./results/DCIS_Precision_Panel_NKI_OncoPrint/oncoPrint_panel_", unique(eventDataFrame$case_control),".pdf"), height = 5)
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
            column_title = paste0(unique(eventDataFrame$case_control), ' in Panel-Seq (n=', length(patients), ')')
  ))
  dev.off()
}

write.csv(samples, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples_V2.txt')

