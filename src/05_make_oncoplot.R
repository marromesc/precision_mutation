library(collapse)
library(ComplexHeatmap)
library(data.table)

setwd('/mnt/albyn/maria/precision_mutation')


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

# Load sample metadata ----------------------------------------------------

#wes
samples_wes <- as.data.frame(fread('./data/WES/DCIS_Precision_WES_All_Samples.txt'))
samples_wes <- samples_wes[samples_wes$qc_normal == 'pass' & samples_wes$qc_pdcis == 'pass',]
samples_wes <- samples_wes[samples_wes$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA)
                           & samples_wes$surgery_final == 'BCS',]
patients_wes <- unique(samples_wes$patient_id); length(patients_wes)
samples_wes$batch <- ifelse(samples_wes$batch %in% c('SLO1', 'SLO2'), 'SLO', ifelse(samples_wes$batch %in% c('NKI1', 'NKI2', 'NKI3', 'NKI4'), 'NKI', ifelse(samples_wes$batch == 'Duke1', 'DUK', ifelse(samples_wes$batch == 'SLO3', 'Melbourne', NA))))

#nki panel
samples_nki_panel <- as.data.frame(fread('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples.txt'))
samples_nki_panel <- samples_nki_panel[samples_nki_panel$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA)
                                       & samples_nki_panel$surgery_final == 'BCS',]
patients_nki_panel <- unique(samples_nki_panel$patient_id); length(patients_nki_panel)
patients_nki_panel <- patients_nki_panel[!(patients_nki_panel %in% patients_wes)]
samples_nki_panel <- samples_nki_panel[samples_nki_panel$patient_id %in% patients_nki_panel,]
samples_nki_panel$batch <- 'NKI'

#kcl panel
samples_kcl_panel <- as.data.frame(fread('./data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_Panel_Sloane_Samples.txt'))
samples_kcl_panel <- samples_kcl_panel[samples_kcl_panel$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA)
                                       & samples_nki_panel$surgery_final == 'BCS',]
patients_kcl_panel <- unique(samples_kcl_panel$patient_id); length(patients_kcl_panel)
patients_kcl_panel <- patients_kcl_panel[!(patients_kcl_panel %in% patients_wes)]
samples_kcl_panel <- samples_kcl_panel[samples_kcl_panel$patient_id %in% patients_kcl_panel,]
samples_kcl_panel$batch <- 'SLO'

patients_all <- unique(c(patients_wes, patients_nki_panel, patients_kcl_panel))
samples <- rbind(samples_wes[,colnames(samples_nki_panel)], samples_nki_panel, samples_kcl_panel[,colnames(samples_nki_panel)])
rownames(samples) <- samples$patient_id

# Load filtered mutations -------------------------------------------------

#wes
eventDataFrame_indel <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds')
eventDataFrame_mutect <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds')

eventDataFrame_mutect <- eventDataFrame_mutect[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')]
eventDataFrame_indel <- eventDataFrame_indel[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                                'er', 'her2', 'grade', 'radiotherapy', 'batch')]

eventDataFrame_all <- rbind(eventDataFrame_mutect, eventDataFrame_indel)
eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_wes,]

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_wes <- eventDataFrame_all

#kcl panel
eventDataFrame_mutect <- readRDS('./data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect_Filtered.rds')

eventDataFrame_all <- eventDataFrame_mutect[,c('ExonicFunc.refGene', 'Chr', 'Start', 'End', 'Ref', 'Alt',
                                               'AF_PDCIS', 'Gene.refGene', 'patient_id', 'first_subseq_event',
                                               'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_kcl_panel & !(eventDataFrame_all$patient_id %in% patients_wes),]
eventDataFrame_all$batch <- 'SLO'

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_kcl_panel,]

eventDataFrame_kcl_panel <- eventDataFrame_all

#nki panel
eventDataFrame_mutect <- readRDS('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.rds')

eventDataFrame_all <- eventDataFrame_mutect[,c('Consequence', 'CHROM', 'POS', 'POS', 'REF', 'ALT',
                                               'vaf_DCIS_DNA', 'Gene.refGene', 'patient_id', 'first_subseq_event',
                                               'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_nki_panel & !(eventDataFrame_all$patient_id %in% patients_wes),]
eventDataFrame_all$batch <- 'NKI'

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_nki_panel,]

eventDataFrame_nki_panel <- eventDataFrame_all

#rbind
eventDataFrame <- rbind(eventDataFrame_kcl_panel, eventDataFrame_nki_panel, eventDataFrame_wes)
eventDataFrame_panel <- rbind(eventDataFrame_nki_panel, eventDataFrame_kcl_panel)

# define cases and controls
eventDataFrame_case <- eventDataFrame[eventDataFrame$case_control == 'case',]
eventDataFrame_control <- eventDataFrame[eventDataFrame$case_control == 'control',]


# Mutation frequency of panel ------------------------------------------------------

GenesPanel <- readRDS('./data/GenesPanel.RDS')
genes_sorted <- readRDS('./data/GenesPanel_sorted.RDS')


# Make oncoprint for wes genes ----------------------------------------------------------

for (status in c('case', 'control')){
  if (status == 'case'){
    eventDataFrame <- eventDataFrame_case
  } else if (status == 'control') {
    eventDataFrame <- eventDataFrame_control
  }
  
  # extract patients
  patients <- samples$patient_id[samples$case_control == status]
  
  # prepare matrix
  geneMatrix <- matrix("",nrow=length(patients),ncol=length(GenesPanel))
  rownames(geneMatrix) <- patients
  colnames(geneMatrix) <- GenesPanel
  ER <- rep("",length(patients));  Her2 <- rep("",length(patients));  Grade <- rep("",length(patients)); RT = rep('', length(patients)); Batch = rep('', length(patients))
  
  for (i in 1:nrow(geneMatrix)) {
    Ind <- which(eventDataFrame$patient_id==rownames(geneMatrix)[i])
    pat_id <- rownames(geneMatrix)[i]
    ER[i] <- samples$er[which(samples$patient_id == pat_id)]
    Her2[i] <- samples$her2[which(samples$patient_id == pat_id)]
    Grade[i] <- samples$grade[which(samples$patient_id == pat_id)]
    RT[i] <- samples$radiotherapy[which(samples$patient_id == pat_id)]
    Batch[i] <-  samples$batch[which(samples$patient_id == pat_id)]
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
  pdf(paste0("./results/oncoPrint_all_", unique(eventDataFrame$case_control),".pdf"), height = 5)
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

write.csv(samples, './data/DCIS_Precision_All.txt')
