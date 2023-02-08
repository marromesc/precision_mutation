# Maria Roman Escorza - 2023 01 18

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

library(collapse)
library(ComplexHeatmap)
library(data.table)
library(dplyr)

source('./lib/oncoPlotDetails.R')
source('./lib/mutCountMatrix.R')

wes_indel_datapath <- './data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered_discovery.rds'
wes_mutect_datapath <- './data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered_discovery.rds'
wes_meta_datapath <- './results/SampleSheet.csv'

census_datapath <- '/mnt/albyn/common/master/cancer_gene_census.csv'
dnds_datapath <- '/mnt/albyn/maria/precision_mutation/results/dmdscv/sel_cv_mutectandpindel.csv'


# Load data ---------------------------------------------------------------

SampleSheet <- as.data.frame(fread(wes_meta_datapath))
SampleSheet <- SampleSheet[SampleSheet$platform=='WES',]

eventDataFrame_mutect <- readRDS(wes_mutect_datapath) %>% dplyr::mutate(Hugo_Symbol=gene.knowngene, Entrez_Gene_Id=patient_id, Center=paste0(batch,'-WES'), NCBI_Build='hg19', Chromosome=chr, Start_Position=start, 
                                                                        End_Position=end, Strand='+', 
                                                                        Variant_Classification=ifelse(exonicfunc.knowngene=='nonsynonymous SNV', 'Missense_Mutation',
                                                                                                      ifelse(exonicfunc.knowngene%in%c('stopgain','stoploss'), 'Nonsense_Mutation',
                                                                                                             ifelse(exonicfunc.knowngene=='.', 'Splice_Site','NoData'))), 
                                                                        Variant_Type='SNP', Reference_Allele=ref_allele, Tumor_Seq_Allele1=alt_allele, Tumor_Seq_Allele2=alt_allele, dbSNP_RS=dbsnp_site,
                                                                        Tumor_Sample_Barcode=sample_name_pdcis, Matched_Norm_Sample_Barcode=sample_name_normal, 
                                                                        t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count, t_vaf=tumor_f, n_vaf=normal_af,
                                                                        SIFT=sift_pred, PolyPhen=polyphen2_hvar_pred, GMAF=x1kg2015aug_max, CLIN_SIG=clinsig, ExAC_AF=exac_all,
                                                                        COSMIC=cosmic70, ESP6500=esp6500siv2_all) %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_Position, 
                                                                                                                                    End_Position, Strand, Variant_Classification, Variant_Type, 
                                                                                                                                    Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, Tumor_Sample_Barcode, 
                                                                                                                                    Matched_Norm_Sample_Barcode, t_depth, t_ref_count, t_alt_count, n_depth, 
                                                                                                                                    n_ref_count, n_alt_count, t_vaf, n_vaf, SIFT, PolyPhen, GMAF, CLIN_SIG, 
                                                                                                                                    ExAC_AF, COSMIC, ESP6500, age_diagnosis, radiotherapy, time_at_risk,
                                                                                                                                    first_subseq_event, case_control, er, pr, her2, grade, rna_id, cnv_id, wes_id)

eventDataFrame_indel <- readRDS(wes_indel_datapath) %>% dplyr::mutate(Hugo_Symbol=gene.knowngene, Entrez_Gene_Id=patient_id, Center=paste0(batch,'-WES'), NCBI_Build='hg19', Chromosome=chr, Start_Position=start, 
                                                                      End_Position=end, Strand='+', 
                                                                      Variant_Classification=ifelse(exonicfunc.knowngene=='frameshift deletion', 'Frame_Shift_Del',
                                                                                                    ifelse(exonicfunc.knowngene%in%c('stopgain','stoploss'), 'Nonsense_Mutation',
                                                                                                           ifelse(exonicfunc.knowngene=='.', 'Splice_Site',
                                                                                                                  ifelse(exonicfunc.knowngene=='nonframeshift deletion', 'In_Frame_Del',
                                                                                                                         ifelse(exonicfunc.knowngene=='nonframeshift insertion', 'In_Frame_Ins',
                                                                                                                                ifelse(exonicfunc.knowngene=='frameshift insertion', 'Frame_Shift_Ins', 'NoData')))))), 
                                                                      Variant_Type=ifelse(ref_allele=='0' & alt_allele=='-', 'DEL',
                                                                                          ifelse(ref_allele=='-' & alt_allele!='-', 'INS', 'NoData')), Reference_Allele=ref_allele, Tumor_Seq_Allele1=alt_allele, Tumor_Seq_Allele2=alt_allele, dbSNP_RS='NoData',
                                                                      Tumor_Sample_Barcode=sample_name_pdcis, Matched_Norm_Sample_Barcode=sample_name_normal, 
                                                                      t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count, t_vaf=tumor_f, n_vaf=normal_af,
                                                                      SIFT=sift_pred, PolyPhen=polyphen2_hvar_pred, GMAF=x1kg2015aug_max, CLIN_SIG=clinsig, ExAC_AF=exac_all,
                                                                      COSMIC=cosmic70, ESP6500=esp6500siv2_all) %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_Position, 
                                                                                                                                  End_Position, Strand, Variant_Classification, Variant_Type, 
                                                                                                                                  Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, Tumor_Sample_Barcode, 
                                                                                                                                  Matched_Norm_Sample_Barcode, t_depth, t_ref_count, t_alt_count, n_depth, 
                                                                                                                                  n_ref_count, n_alt_count, t_vaf, n_vaf, SIFT, PolyPhen, GMAF, CLIN_SIG, 
                                                                                                                                  ExAC_AF, COSMIC, ESP6500, age_diagnosis, radiotherapy, time_at_risk,
                                                                                                                                  first_subseq_event, case_control, er, pr, her2, grade, rna_id, cnv_id, wes_id)

eventDataFrame <- rbind(eventDataFrame_mutect, eventDataFrame_indel)

eventDataFrame <- eventDataFrame[!is.na(eventDataFrame$case_control),]
eventDataFrame <- eventDataFrame[eventDataFrame$Entrez_Gene_Id %in% SampleSheet$patient_id,]

census <- read_csv(census_datapath)[[1]]
dnds <- read_csv(dnds_datapath)
dnds <- dnds[dnds$qglobal_cv<=0.05,][[1]]

eventDataFrame <- eventDataFrame[eventDataFrame$Hugo_Symbol%in%c(census,dnds),]

eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c("Splice_Site")] <- "splicing"
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c('Missense_Mutation')] <- "missense"
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c("In_Frame_Del","In_Frame_Ins")] <- "inframe_indel"
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins")] <- "frameshift" 
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c("Nonsense_Mutation")] <- "nonsense"

GenesPanel_ca <- table(eventDataFrame$Hugo_Symbol[eventDataFrame$case_control=='case'])
GenesPanel_ca <- GenesPanel_ca[GenesPanel_ca>6]
GenesPanel_co <- table(eventDataFrame$Hugo_Symbol[eventDataFrame$case_control=='control'])
GenesPanel_co <- GenesPanel_co[GenesPanel_co>6]
GenesPanel <- unique(c(names(GenesPanel_ca), names(GenesPanel_co)))


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
    Ind <- which(eventDataFrame$Entrez_Gene_Id==rownames(geneMatrix)[i])
    pat_id <- rownames(geneMatrix)[i]
    ER[i] <- SampleSheet$er[which(SampleSheet$patient_id == pat_id)]
    Her2[i] <- SampleSheet$her2[which(SampleSheet$patient_id == pat_id)]
    Grade[i] <- SampleSheet$grade[which(SampleSheet$patient_id == pat_id)]
    RT[i] <- SampleSheet$radiotherapy[which(SampleSheet$patient_id == pat_id)]
    Batch[i] <-  SampleSheet$cohort[which(SampleSheet$patient_id == pat_id)]
    for (j in 1:ncol(geneMatrix)) {
      Ind <- which(eventDataFrame$Entrez_Gene_Id==rownames(geneMatrix)[i] & eventDataFrame$Hugo_Symbol==colnames(geneMatrix)[j])
      if (length(Ind)==1) {
        geneMatrix[i,j] <- eventDataFrame$Consequence[Ind]
      }
      if (length(Ind)>1) {
        if (length(unique(eventDataFrame$Entrez_Gene_Id[Ind]))==1) {  
          geneMatrix[i,j] <- "multi_hit"
        }
        else {
          geneMatrix[i,j] <- eventDataFrame$Consequence[Ind[1]]
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
  pdf(paste0("./results/oncoPrint_wes_", status,".pdf"), height = 6.5)
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
                  column_title = paste0(status, ' (n=', length(patients), ')')
  ))
  dev.off()
}
