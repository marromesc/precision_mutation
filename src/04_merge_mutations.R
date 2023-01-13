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
eventDataFrame_indel <- readRDS(wes_indel_datapath)
eventDataFrame_mutect <- readRDS(wes_mutect_datapath)

eventDataFrame_mutect <- eventDataFrame_mutect[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')]
eventDataFrame_indel <- eventDataFrame_indel[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                                'er', 'her2', 'grade', 'radiotherapy', 'batch')]

eventDataFrame_all <- rbind(eventDataFrame_mutect, eventDataFrame_indel)
eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% SampleSheet$patient_id[SampleSheet$platform == 'WES'],]

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_wes <- eventDataFrame_all

#kcl panel
eventDataFrame_mutect <- readRDS(kcl_datapath)

eventDataFrame_all <- eventDataFrame_mutect[,c('ExonicFunc.refGene', 'Chr', 'Start', 'End', 'Ref', 'Alt',
                                               'AF_PDCIS', 'Gene.refGene', 'patient_id', 'first_subseq_event',
                                               'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% SampleSheet$patient_id[SampleSheet$platform == 'PanelSeq'],]
eventDataFrame_all$batch <- 'SLO'

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_kcl_panel <- eventDataFrame_all

#nki panel
eventDataFrame_mutect <- readRDS(nki_datapath)

eventDataFrame_all <- eventDataFrame_mutect[,c('Consequence', 'CHROM', 'POS', 'POS', 'REF', 'ALT',
                                               'vaf_DCIS_DNA', 'Gene.refGene', 'patient_id', 'first_subseq_event',
                                               'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% SampleSheet$patient_id[SampleSheet$platform == 'PanelSeq'],]
eventDataFrame_all$batch <- 'NKI'

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_nki_panel <- eventDataFrame_all


# Merge platforms ---------------------------------------------------------

#rbind
eventDataFrame <- rbind(eventDataFrame_kcl_panel, eventDataFrame_nki_panel, eventDataFrame_wes)


# Make maf file -----------------------------------------------------------

SampleSheet_maf <- SampleSheet %>% mutate(Tumor_Sample_Barcode = patient_id, Grade = grade, days_to_last_followup = fup_months, Overall_Survival_Status = ifelse(case_control == 'case', 1, 0) )
eventDataFrame$NCBI_Build <- 37
eventDataFrame$Strand <- '+'
eventDataFrame$Variant_Classification <- ifelse(eventDataFrame$exonicfunc.knowngene == '.', 'Splice_Site', 
                                                ifelse(eventDataFrame$exonicfunc.knowngene %in% c('nonsynonymous SNV', 'missense_variant'), 'Missense_Mutation',
                                                       ifelse(eventDataFrame$exonicfunc.knowngene == 'frameshift deletion', 'Frame_Shift_Del',
                                                              ifelse(eventDataFrame$exonicfunc.knowngene %in% c('nonframeshift substitution', 'nonframeshift deletion'), 'In_Frame_Del',
                                                                     ifelse(eventDataFrame$exonicfunc.knowngene %in% c('stopgain'), 'Nonsense_Mutation',
                                                                            ifelse(eventDataFrame$exonicfunc.knowngene == 'frameshift insertion', 'Frame_Shift_Ins',
                                                                                   ifelse(eventDataFrame$exonicfunc.knowngene == 'nonframeshift insertion', 'In_Frame_Ins', NA)))))))

eventDataFrame$Variant_Type <- ifelse(eventDataFrame$Variant_Classification %in% c('Frame_Shift_Ins', 'In_Frame_Ins'), 'INS',
                                      ifelse(eventDataFrame$Variant_Classification %in% c('Frame_Shift_Del', 'In_Frame_Del'), 'INS', 'SNP'))
maf <- eventDataFrame %>% dplyr::select(Hugo_Symbol = gene.knowngene, Entrez_Gene_Id = patient_id, Center = batch, NCBI_Build, Chromosome = chr, Start_Position = start, End_Position = end, Variant_Classification, Variant_Type, Reference_Allele = ref_allele, Tumor_Seq_Allele1 = ref_allele, Tumor_Seq_Allele2 = alt_allele, Tumor_Sample_Barcode = patient_id)

write.table(SampleSheet_maf, './results/SampleSheet_maf.tsv', row.names = F, quote = F, sep='\t')
write.table(maf, './results/Filtered_Mutations_Compiled.maf', row.names = F, quote = F, sep='\t')


# Rename mutations and export ---------------------------------------------
# filter mutations in panel genes
eventDataFrame <- eventDataFrame[eventDataFrame$gene.knowngene %in% GenesPanel,]

# rename mutations
eventDataFrame$exonicfunc.knowngene[eventDataFrame$exonicfunc.knowngene %in% c(".")] <- "splicing"
eventDataFrame$exonicfunc.knowngene[eventDataFrame$exonicfunc.knowngene %in% c("nonsynonymous SNV", 'missense_variant')] <- "missense"
eventDataFrame$exonicfunc.knowngene[eventDataFrame$exonicfunc.knowngene %in% c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution")] <- "inframe_indel"
eventDataFrame$exonicfunc.knowngene[eventDataFrame$exonicfunc.knowngene %in% c("frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift" 
eventDataFrame$exonicfunc.knowngene[eventDataFrame$exonicfunc.knowngene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"

write.csv(eventDataFrame, './results/Filtered_Mutations_Compiled.csv', row.names = F)
