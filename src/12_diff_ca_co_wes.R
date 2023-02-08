# Maria Roman Escorza - 2023 01 18

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

system('mkdir ./results/diff_mut_wes')

library(readr)
library(data.table)
library(dplyr)

source('./lib/forestPlotFisher.R')
source('./lib/diffMut.R')
source('./lib/mutCountMatrix.R')
source('./lib/addCN2Muts.R')

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

GenesPanel_ca <- table(eventDataFrame$Hugo_Symbol[eventDataFrame$case_control=='case'])
GenesPanel_ca <- GenesPanel_ca[GenesPanel_ca>6]
GenesPanel_co <- table(eventDataFrame$Hugo_Symbol[eventDataFrame$case_control=='control'])
GenesPanel_co <- GenesPanel_co[GenesPanel_co>6]
GenesPanel <- unique(c(names(GenesPanel_ca), names(GenesPanel_co)))


# Differential analysis ---------------------------------------------------

geneMatrix <- mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, GenesPanel)
fisherTable <- diffMut(geneMatrix, SampleSheet)
write.table(fisherTable, './results/diff_mut_wes/diff_mut_overall.csv', sep = ',', row.names = F)

#forest plot
if (length(which(fisherTable$adjPval < 0.05))>0){
  pdf('./results/diff_mut_wes/diff_mut_overall.pdf')
  forestPlotFisher(fisherTable, SampleSheet, fdr=0.05)
  dev.off()
}


# Differential analysis RT ---------------------------------------------------

for (rt in c(0, 1)){
  geneMatrix_i <- geneMatrix[SampleSheet$patient_id[SampleSheet$radiotherapy == rt],]
  fisherTable <- diffMut(geneMatrix_i, SampleSheet[SampleSheet$radiotherapy == rt,])
  write.table(fisherTable, paste0('./results/diff_mut_wes/diff_mut_RT_', rt, '.csv'), sep = ',', row.names = F)
  
  if (length(which(fisherTable$adjPval < 0.05))>0){
    pdf( paste0('./results/diff_mut_wes/diff_mut_RT_', rt, '.pdf'))
    forestPlotFisher(fisherTable, SampleSheet = SampleSheet[SampleSheet$radiotherapy == rt,], fdr=0.05)
    dev.off()
  }
}


# Differential analysis Her2 ---------------------------------------------------

for (her2 in c(0, 1)){
  geneMatrix_i <- geneMatrix[SampleSheet$patient_id[SampleSheet$her2 == her2],]
  fisherTable <- diffMut(geneMatrix_i, SampleSheet[SampleSheet$her2 == her2,])  
  write.table(fisherTable, paste0('./results/diff_mut_wes/diff_mut_HER2_', her2, '.csv'), sep = ',', row.names = F)
  
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
  write.table(fisherTable, paste0('./results/diff_mut_wes/diff_mut_ER_', er, '.csv'), sep = ',', row.names = F)
  
  if (length(which(fisherTable$adjPval < 0.05))>0){
    pdf( paste0('./results/diff_mut_wes/diff_mut_ER_', er, '.pdf'))
    forestPlotFisher(fisherTable, SampleSheet = SampleSheet[SampleSheet$er == er,], fdr=0.05)
    dev.off()
  }
}
