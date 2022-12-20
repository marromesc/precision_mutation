library(data.table)
library(readr)
library(ggplot2)
library(ggpubr)

setwd('/mnt/albyn/maria/precision_mutation')

# Load sample metadata ----------------------------------------------------

#wes
samples_wes <- as.data.frame(fread('./data/WES/DCIS_Precision_WES_All_Samples.txt'))
samples_wes <- samples_wes[samples_wes$qc_normal == 'pass' & samples_wes$qc_pdcis == 'pass',]
samples_wes <- samples_wes[samples_wes$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA),]
patients_wes <- unique(samples_wes$patient_id); length(patients_wes)

#nki panel
samples_nki_panel <- as.data.frame(fread('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples.txt'))
samples_nki_panel <- samples_nki_panel[samples_nki_panel$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA),]
patients_nki_panel <- unique(samples_nki_panel$patient_id); length(patients_nki_panel)

#kcl panel
samples_kcl_panel <- as.data.frame(fread('./data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_Panel_Sloane_Samples.txt'))
samples_kcl_panel <- samples_kcl_panel[samples_kcl_panel$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA),]
patients_kcl_panel <- unique(samples_kcl_panel$patient_id); length(patients_kcl_panel)

patients_all <- unique(c(patients_wes, patients_nki_panel, patients_kcl_panel))
  
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
eventDataFrame_all$batch <- 'SLO_Panel'

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_kcl_panel <- eventDataFrame_all

#nki panel
eventDataFrame_mutect <- readRDS('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.rds')

eventDataFrame_all <- eventDataFrame_mutect[,c('Consequence', 'CHROM', 'POS', 'POS', 'REF', 'ALT',
                                               'vaf_DCIS_DNA', 'Gene.refGene', 'patient_id', 'first_subseq_event',
                                               'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_nki_panel & !(eventDataFrame_all$patient_id %in% patients_wes),]
eventDataFrame_all$batch <- 'NKI_Panel'

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_nki_panel <- eventDataFrame_all

#rbind
eventDataFrame_all <- rbind(eventDataFrame_kcl_panel, eventDataFrame_nki_panel, eventDataFrame_wes)


# Mutation count matrix ---------------------------------------------------

GenesPanel <- readRDS('./data/GenesPanel.RDS')

geneMatrix <- matrix("",nrow=length(patients_all),ncol=length(GenesPanel))
rownames(geneMatrix) <- patients_all
colnames(geneMatrix) <- GenesPanel

for (i in 1:nrow(geneMatrix)) {
  Ind <- which(eventDataFrame_all$patient_id==rownames(geneMatrix)[i])
  pat_id <- rownames(geneMatrix)[i]
  for (j in 1:ncol(geneMatrix)) {
    Ind <- which(eventDataFrame_all$patient_id==rownames(geneMatrix)[i] & eventDataFrame_all$gene.knowngene==colnames(geneMatrix)[j])
    if (length(Ind)==1) {
      geneMatrix[i,j] <- eventDataFrame_all$exonicfunc.knowngene[Ind]
    }
    if (length(Ind)>1) {
      if (length(unique(eventDataFrame_all$patient_id[Ind]))==1) {  
        geneMatrix[i,j] <- "multi_hit"
      }
      else {
        geneMatrix[i,j] <- eventDataFrame_all$exonicfunc.knowngene[Ind[1]]
      }
    }
  }
}

geneMatrix <- ifelse(geneMatrix == "", 0, 1)
no_mutated <- as.data.frame(rowSums(geneMatrix))
#no_mutated <- no_mutated[no_mutated==0]


# Load cellularity  -------------------------------------------------------

openclinica <- as.data.frame(readxl::read_xlsx('./data/updated_20221114_clinicaldata_caco_genomic_samples.xlsx'))
openclinica2 <- as.data.frame(readxl::read_xlsx('./data/updated_20221114_clinicaldata_caco_genomic_samples.xlsx', sheet = 2))
openclinica2 <- openclinica2[openclinica2$redcap_event == 'PRI',]
openclinica <- merge(openclinica, openclinica2)
openclinica <- openclinica[openclinica$patient_id %in% eventDataFrame_all$patient_id,]
openclinica <- openclinica[!is.na(openclinica$cnv_id),]

cellularity <- as.data.frame(read_tsv('/mnt/albyn/maria/copynumber/precision_copynumber/results/segmentation/2N/fitpicker_2N.tsv'))
cellularity <- cellularity[cellularity$sample %in% openclinica$cnv_id,c('sample', 'likely_fit')]
cellularity <- merge(cellularity, openclinica, by.x='sample', by.y='cnv_id'); dim(cellularity)
cellularity <- merge(cellularity, no_mutated, by.x='patient_id', by.y=0)
colnames(cellularity)[ncol(cellularity)] <- 'n_muts'
cellularity$mutations <- ifelse(cellularity$n_muts == 0, 'no_mutations', 'mutations')

corResult <- cor.test(cellularity$likely_fit, cellularity$n_muts, method = 'pearson')

mainTitle <- paste0("Correlation = ",
                    round(corResult[["estimate"]], 2),
                    "   P-value = ",
                    signif(corResult[["p.value"]], 2))

pdf('./results/correlation_cellularity_n_muts.pdf')
ggplot(cellularity, aes(x=likely_fit, y=n_muts)) + 
  geom_point()+
  geom_smooth(method=lm) +
  labs(title=mainTitle)

ggboxplot(cellularity, x = "mutations", y = "likely_fit", palette = "jco",
               add = "jitter") + stat_compare_means() 
dev.off()



# CII score ---------------------------------------------------------------

cii <- as.data.frame(read_tsv('/mnt/albyn/maria/precision-CaseControl/Tables/cii.tsv'))
cii <- cii[cii$id %in% cellularity$sample,]
cellularity <- merge(cellularity, cii, by.x='sample', by.y='id')

corResult <- cor.test(cellularity$cii, cellularity$n_muts, method = 'pearson')

mainTitle <- paste0("Correlation = ",
                    round(corResult[["estimate"]], 2),
                    "   P-value = ",
                    signif(corResult[["p.value"]], 2))

pdf('./results/correlation_cii_n_muts.pdf')
ggplot(cellularity, aes(x=cii, y=n_muts)) + 
  geom_point()+
  geom_smooth(method=lm) +
  labs(title=mainTitle)

ggboxplot(cellularity, x = "mutations", y = "cii", palette = "jco",
          add = "jitter") + stat_compare_means()

ggboxplot(cellularity, x = "n_muts", y = "cii", palette = "jco",
          add = "jitter") + stat_compare_means()
dev.off()


