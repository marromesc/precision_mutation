# Maria Roman Escorza - 2022 12 

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')
system('mkdir ./results/tmb')

library(maftools)
library(ggplot2)
library(ggpubr)
library(readr)

wes_indel_datapath <- './data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds'
wes_mutect_datapath <- './data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds'
wes_meta_datapath <- './results/SampleSheet.csv'


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

eventDataFrame_wes <- rbind(eventDataFrame_mutect, eventDataFrame_indel)

eventDataFrame_wes <- eventDataFrame_wes[!is.na(eventDataFrame_wes$case_control),]
eventDataFrame_wes <- eventDataFrame_wes[eventDataFrame_wes$Entrez_Gene_Id %in% SampleSheet$patient_id,]

maf_case <- read.maf(eventDataFrame_wes[eventDataFrame_wes$Entrez_Gene_Id %in% SampleSheet$patient_id[SampleSheet$case_control=='case'],])
maf_control <- read.maf(eventDataFrame_wes[eventDataFrame_wes$Entrez_Gene_Id %in% SampleSheet$patient_id[SampleSheet$case_control=='control'],])


# Compute TMB -------------------------------------------------------------

pdf('./results/tmb/tmb_case.pdf')
tmb_case <- tmb(maf_case)
dev.off()

pdf('./results/tmb/tmb_control.pdf')
tmb_control <- tmb(maf_control)
dev.off()

tmb_case$case_control <- 'case'
tmb_control$case_control <- 'control'
tmb <- rbind(tmb_case, tmb_control)

pdf('./results/tmb/total_tmb.pdf', width = 4, height = 4, colormodel = "cmyk", useDingbats = TRUE)
ggplot(tmb, aes(x=case_control, y=total, color=case_control)) +
  geom_boxplot() +
  geom_point(aes(color = case_control), position = position_jitterdodge()) +
  stat_compare_means() +
  theme(legend.position="none", plot.title = element_text(lineheight=3, face="bold", color="black", size=15)) +
  scale_color_manual(values=c("#648FFF", '#DA0E91')) +
  labs(y = "Total TMB", x="") +
  theme(axis.title.y = element_text(size = rel(1), angle = 90), axis.line = element_line(colour = 'black'), panel.background = element_blank()) +
  scale_x_discrete(labels = c('control' = 'Control', 'case' = 'INV-Case')) +
  stat_summary(fun = mean, geom='point', fill = 'black', shape = 18, size = 3, color='black')
dev.off() 

pdf('./results/tmb/total_perMB.pdf', width = 4, height = 4, colormodel = "cmyk", useDingbats = TRUE)
ggplot(tmb, aes(x=case_control, y=total_perMB, color=case_control)) +
  geom_boxplot() +
  geom_point(aes(color = case_control), position = position_jitterdodge()) +
  stat_compare_means() +
  theme(legend.position="none", plot.title = element_text(lineheight=3, face="bold", color="black", size=15)) +
  scale_color_manual(values=c("#648FFF", '#DA0E91')) +
  labs(y = "Total TMB per MB (log)", x="") +
  theme(axis.title.y = element_text(size = rel(1), angle = 90), axis.line = element_line(colour = 'black'), panel.background = element_blank()) +
  scale_x_discrete(labels = c('control' = 'Control', 'case' = 'INV-Case')) +
  stat_summary(fun = mean, geom='point', fill = 'black', shape = 18, size = 3, color='black')
dev.off() 

