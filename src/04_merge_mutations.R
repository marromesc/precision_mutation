# Maria Roman Escorza - 2023 01 04

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

library(data.table)
library(dplyr)

wes_meta_datapath <- './data/WES/DCIS_Precision_WES_All_Samples.txt'
nki_panel_meta_datapath <- './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples.txt'
slo_panel_meta_datapath <- './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_Panel_Sloane_Samples.txt'

wes_indel_datapath <- './data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds'
wes_mutect_datapath <- './data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds'
kcl_datapath <- './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect_Filtered.rds'
nki_datapath <- './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.rds'

GenesPanel <- readRDS('./data/GenesPanel.RDS')


# Load sample metadata ----------------------------------------------------

#wes
SampleSheet_wes <- as.data.frame(fread(wes_meta_datapath))
SampleSheet_wes <- SampleSheet_wes[SampleSheet_wes$qc_normal == 'pass' & SampleSheet_wes$qc_pdcis == 'pass',]
SampleSheet_wes$cohort <- ifelse(SampleSheet_wes$batch %in% c('SLO1', 'SLO2'), 'SLO', ifelse(SampleSheet_wes$batch %in% c('NKI1', 'NKI2', 'NKI3', 'NKI4'), 'NKI', ifelse(SampleSheet_wes$batch == 'Duke1', 'DUK', ifelse(SampleSheet_wes$batch == 'SLO3', 'Melbourne', NA))))
SampleSheet_wes$platform <- 'WES'

#nki panel
SampleSheet_nki_panel <- as.data.frame(fread(nki_panel_meta_datapath))
SampleSheet_nki_panel <- SampleSheet_nki_panel[is.na(SampleSheet_nki_panel$wes_id),]
SampleSheet_nki_panel$batch <- 'NKI1'
SampleSheet_nki_panel$cohort <- 'NKI'
SampleSheet_nki_panel$platform <- 'PanelSeq'

#kcl panel
SampleSheet_slo_panel <- as.data.frame(fread(slo_panel_meta_datapath))
SampleSheet_slo_panel <- SampleSheet_slo_panel[is.na(SampleSheet_slo_panel$wes_id),]
SampleSheet_slo_panel$batch <- 'SLO1'
SampleSheet_slo_panel$cohort <- 'SLO'
SampleSheet_slo_panel$platform <- 'PanelSeq'

#merge metadata from different sequencing platforms
SampleSheet <- rbind(SampleSheet_wes[,colnames(SampleSheet_nki_panel)], SampleSheet_nki_panel, SampleSheet_slo_panel[,colnames(SampleSheet_nki_panel)])
rownames(SampleSheet) <- SampleSheet$patient_id
SampleSheet <- SampleSheet[SampleSheet$surgery_final == 'BCS' & SampleSheet$first_subseq_event %in% c('ipsilateral IBC', NA, 'death'),]
write.csv(SampleSheet, './results/SampleSheet.csv', row.names = F)


# Load filtered mutations -------------------------------------------------

#wes
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

#kcl panel
eventDataFrame_kcl_panel <- readRDS(kcl_datapath) %>% dplyr::mutate(Hugo_Symbol=Gene.refGene, Entrez_Gene_Id=patient_id, Center='SLO-Panel', NCBI_Build='hg19', Chromosome=Chr, Start_Position=Start, 
                                                                    End_Position=End, Strand='+', 
                                                                    Variant_Classification=ifelse(ExonicFunc.refGene=='nonsynonymous SNV', 'Missense_Mutation',
                                                                                                  ifelse(ExonicFunc.refGene%in%c('stopgain','stoploss'), 'Nonsense_Mutation',
                                                                                                         ifelse(ExonicFunc.refGene=='.', 'Splice_Site',
                                                                                                                ifelse(ExonicFunc.refGene=='frameshift deletion', 'Frame_Shift_Del',
                                                                                                                       ifelse(ExonicFunc.refGene=='frameshift insertion', 'Frame_Shift_Ins',
                                                                                                                              ifelse(ExonicFunc.refGene=='nonframeshift deletion', 'In_Frame_Del',
                                                                                                                                     ifelse(ExonicFunc.refGene%in%c('nonframeshift insertion','nonframeshift substitution'), 'In_Frame_Ins', 'NoData'))))))), 
                                                                    Variant_Type=ifelse(nchar(Ref)==2 & nchar(Alt)==2, 'DNP',
                                                                                        ifelse(Alt=='-', 'DEL',
                                                                                               ifelse(Ref=='-', 'INS', 'SNP'))), Reference_Allele=Ref, Tumor_Seq_Allele1=Alt, Tumor_Seq_Allele2=Alt, dbSNP_RS='NoData',
                                                                    Tumor_Sample_Barcode='NoData', Matched_Norm_Sample_Barcode='NoData', 
                                                                    t_depth='NoData', t_ref_count='NoData', t_alt_count='NoData', n_depth='NoData', n_ref_count='NoData', n_alt_count='NoData', t_vaf=AF_PDCIS, n_vaf=AF_NOR,
                                                                    SIFT=SIFT_pred, PolyPhen=Polyphen2_HVAR_pred, GMAF=ALL.sites.2015_08, CLIN_SIG=CLINSIG, ExAC_AF=ExAC_ALL,
                                                                    COSMIC=cosmic70, ESP6500=esp6500siv2_all) %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_Position, 
                                                                                                       End_Position, Strand, Variant_Classification, Variant_Type, 
                                                                                                       Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, Tumor_Sample_Barcode, 
                                                                                                       Matched_Norm_Sample_Barcode, t_depth, t_ref_count, t_alt_count, n_depth, 
                                                                                                       n_ref_count, n_alt_count, t_vaf, n_vaf, SIFT, PolyPhen, GMAF, CLIN_SIG, 
                                                                                                       ExAC_AF, COSMIC, ESP6500, age_diagnosis, radiotherapy, time_at_risk,
                                                                                                       first_subseq_event, case_control, er, pr, her2, grade, rna_id, cnv_id, wes_id)

#nki panel
eventDataFrame_nki_panel <- readRDS(nki_datapath) %>% dplyr::mutate(Hugo_Symbol=SYMBOL, Entrez_Gene_Id=patient_id, Center='NKI-Panel', NCBI_Build='hg19', Chromosome=CHROM, Start_Position=POS, 
                                                                    End_Position=POS, Strand='+', 
                                                                    Variant_Classification=ifelse(ExonicFunc.refGene=='nonsynonymous SNV', 'Missense_Mutation',
                                                                                                  ifelse(ExonicFunc.refGene%in%c('stopgain','stoploss'), 'Nonsense_Mutation',
                                                                                                          ifelse(ExonicFunc.refGene%in%c('.',''), 'Splice_Site',
                                                                                                                 ifelse(ExonicFunc.refGene=='frameshift deletion', 'Frame_Shift_Del',
                                                                                                                        ifelse(ExonicFunc.refGene=='frameshift insertion', 'Frame_Shift_Ins',
                                                                                                                               ifelse(ExonicFunc.refGene=='nonframeshift deletion', 'In_Frame_Del',
                                                                                                                                     ifelse(ExonicFunc.refGene%in%c('nonframeshift insertion','nonframeshift substitution'), 'In_Frame_Ins', 'NoData'))))))), 
                                                                    Variant_Type=ifelse(nchar(REF)==2 & nchar(ALT)==2, 'DNP',
                                                                                        ifelse(ALT=='-', 'DEL',
                                                                                               ifelse(REF=='-', 'INS', 'SNP'))), Reference_Allele=REF, Tumor_Seq_Allele1=ALT, Tumor_Seq_Allele2=ALT, dbSNP_RS='NoData',
                                                                    Tumor_Sample_Barcode='NoData', Matched_Norm_Sample_Barcode='NoData', 
                                                                    t_depth='NoData', t_ref_count='NoData', t_alt_count='NoData', n_depth='NoData', n_ref_count='NoData', n_alt_count='NoData', t_vaf=as.numeric(vaf_DCIS_DNA)/100, n_vaf='NoData',
                                                                    SIFT=SIFT_pred, PolyPhen=Polyphen2_HVAR_pred, GMAF=ALL.sites.2015_08, CLIN_SIG=CLIN_SIG, ExAC_AF=ExAC_ALL,
                                                                    COSMIC=cosmic70, ESP6500=esp6500siv2_all) %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_Position, 
                                                                                                                                End_Position, Strand, Variant_Classification, Variant_Type, 
                                                                                                                                Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, Tumor_Sample_Barcode, 
                                                                                                                                Matched_Norm_Sample_Barcode, t_depth, t_ref_count, t_alt_count, n_depth, 
                                                                                                                                n_ref_count, n_alt_count, t_vaf, n_vaf, SIFT, PolyPhen, GMAF, CLIN_SIG, 
                                                                                                                                ExAC_AF, COSMIC, ESP6500, age_diagnosis, radiotherapy, time_at_risk,
                                                                                                                                first_subseq_event, case_control, er, pr, her2, grade, rna_id, cnv_id, wes_id)


# Merge platforms ---------------------------------------------------------

eventDataFrame <- rbind(eventDataFrame_kcl_panel, eventDataFrame_nki_panel, eventDataFrame_wes)
eventDataFrame <- eventDataFrame[eventDataFrame$Entrez_Gene_Id %in% SampleSheet$patient_id,]

# Rename mutations and export ---------------------------------------------
# filter mutations in panel genes
eventDataFrame <- eventDataFrame[eventDataFrame$Hugo_Symbol %in% GenesPanel,]

# rename mutations
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c("Splice_Site")] <- "splicing"
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c('Missense_Mutation')] <- "missense"
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c("In_Frame_Del","In_Frame_Ins")] <- "inframe_indel"
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins")] <- "frameshift" 
eventDataFrame$Consequence[eventDataFrame$Variant_Classification %in% c("Nonsense_Mutation")] <- "nonsense"

write.csv(eventDataFrame, './results/Filtered_Mutations_Compiled.csv', row.names = F)
