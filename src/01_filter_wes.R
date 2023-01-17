# Maria Roman Escorza - 2022 11 30 

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

library(tidyr)
library(dplyr)
library(stringr)
library(gridExtra)

source('./lib/readMetadata.R')
source('./lib/readMutect.R')
source('./lib/MyNumeric.R')
source('./lib/GetCosmicNumber.R')
source('./lib/filterByBed.R')

openclinica_datapath <- './data/updated_20221114_clinicaldata_caco_genomic_samples.xlsx'
qc_datapath <- './data/WES/20221102_DCIS_Precision_WES_sampleInfo_by_XiaogangWu.xlsx'
wes_datapath <- '/mnt/albyn/maria/prj_precision/wes/vcf_annovar/'

bed <- read.delim('./data/WES/Exome-Agilent_V4.bed', sep = "\t", header = FALSE, as.is = TRUE)

# read openclinica master sheet
openclinica <- readMetadata(openclinica_datapath)

# read qc master sheet
sloane_qc <- readxl::read_xlsx(qc_datapath, sheet = 'SLO')
nki_qc <- readxl::read_xlsx(qc_datapath, sheet = 'NKI') ; nki_qc <- nki_qc[!(nki_qc$note %in% c("No sample of primary DCIS, only DCIS recurrence.", "\"not usable\"")),] # remove not usable samples
duke_qc <- readxl::read_xlsx(qc_datapath, sheet = 'Duke')
colnames <- colnames(sloane_qc)[colnames(sloane_qc) %in% colnames(nki_qc)]
wes_qc <- rbind(sloane_qc[,colnames], nki_qc[,colnames], duke_qc[,colnames])


# Coverage Quality Control ------------------------------------------------

# Perform QC (Total reads >= 10M AND median target coverage >= 40 for normal and >= 80 for tumor samples)
wes_qc$qc <- ifelse(wes_qc$`total_reads(million)` >= 10 & wes_qc$tissue_pathology == 'Normal' & wes_qc$median_Target_Cov >= 40, 'pass', 
                ifelse(wes_qc$`total_reads(million)` >= 10 & wes_qc$tissue_pathology != 'Normal' & wes_qc$median_Target_Cov >= 80, 'pass', 'fail'))

#add exception for DCIS00402, ceiling
wes_qc$qc[wes_qc$sampleid == 'DCIS00402Normal-N'] <- 'pass'

# change id for icicle sample
wes_qc$precision_patient_id <- ifelse(wes_qc$precision_patient_id %in% c('IT02006 IR', 'IT02006 N', 'IT02006 P'), 'IT02006', wes_qc$precision_patient_id)
wes_qc$sample_label <- ifelse(wes_qc$precision_patient_id == 'IT02006 P', 'IT02006P', wes_qc$sample_label)

# ignoring failed run due to the issue of sequencing library preparation kits
wes_qc <- wes_qc[!(wes_qc$runid == '211111_A00482_0262_AH3KJWDSX3' & wes_qc$batch == 'SLO3'),]

# removing duke2 dataset
wes_qc <- wes_qc[wes_qc$batch != 'Duke2',]

# add qc data to openclinica sheet
openclinica <- openclinica[openclinica$patient_id %in% wes_qc$precision_patient_id,]; dim(openclinica)
openclinica$batch <- NA

openclinica$sampleid_normal <- NA
openclinica$sample_label_normal <- NA
openclinica$sample_name_normal <- NA
openclinica$sample_runid_normal <- NA
openclinica$total_reads_million_normal <- NA
openclinica$median_cov_normal <- NA
openclinica$mean_cov_normal <- NA
openclinica$qc_normal <- NA

openclinica$sampleid_pdcis <- NA
openclinica$sample_label_pdcis <- NA
openclinica$sample_name_pdcis <- NA
openclinica$sample_runid_pdcis <- NA
openclinica$total_reads_million_pdcis <- NA
openclinica$median_cov_pdcis <- NA
openclinica$mean_cov_pdcis <- NA
openclinica$qc_pdcis <- NA

for (i in 1:nrow(openclinica)){
  pat_id <- openclinica$patient_id[i]
  
  Ind_N <- which(wes_qc$precision_patient_id == pat_id & wes_qc$event_type == 'normal')  #check if there's normal sample for this patient
  
  if (length(Ind_N) != 0){
    openclinica$batch[i] <- wes_qc$batch[Ind_N]
    openclinica$sampleid_normal[i] <- wes_qc$sampleid[Ind_N]
    openclinica$sample_label_normal[i] <- wes_qc$sample_label[Ind_N]
    openclinica$sample_name_normal[i] <- wes_qc$sample_name[Ind_N]
    openclinica$sample_runid_normal[i] <- wes_qc$runid[Ind_N]
    openclinica$total_reads_million_normal[i] <- wes_qc$'total_reads(million)'[Ind_N]
    openclinica$median_cov_normal[i] <- wes_qc$'median_Target_Cov'[Ind_N]
    openclinica$mean_cov_normal[i] <- wes_qc$'Mean_Target_COV'[Ind_N]
    openclinica$qc_normal[i] <- wes_qc$qc[Ind_N]
  }
  
  Ind_P <- which(wes_qc$precision_patient_id == pat_id & wes_qc$event_type == 'primary') #check if there's primary sample for this patient
  
  if (length(Ind_P) != 0){
    openclinica$sampleid_pdcis[i] <- wes_qc$sampleid[Ind_P]
    openclinica$sample_label_pdcis[i] <- wes_qc$sample_label[Ind_P]
    openclinica$sample_name_pdcis[i] <- wes_qc$sample_name[Ind_P]
    openclinica$sample_runid_pdcis[i] <- wes_qc$runid[Ind_P]
    openclinica$total_reads_million_pdcis[i] <- wes_qc$'total_reads(million)'[Ind_P]
    openclinica$median_cov_pdcis[i] <- wes_qc$'median_Target_Cov'[Ind_P]
    openclinica$mean_cov_pdcis[i] <- wes_qc$'Mean_Target_COV'[Ind_P]
    openclinica$qc_pdcis[i] <- wes_qc$qc[Ind_P]
  }
}

# remove patients with no primary or normal
openclinica <- openclinica[!(is.na(openclinica$sampleid_normal) | is.na(openclinica$sampleid_pdcis)),]; dim(openclinica)
write.table(openclinica, './data/WES/DCIS_Precision_WES_All_Samples.txt', sep = '\t', quote = FALSE, row.names = FALSE)

# filter cases and controls
openclinica <- openclinica[openclinica$first_subseq_event %in% c('ipsilateral DCIS', 'ipsilateral IBC', 'NA', 'death') & openclinica$surgery_final == 'BCS',]; dim(openclinica)

# filter samples which passed qc
openclinica <- openclinica[openclinica$qc_normal == 'pass' & openclinica$qc_pdcis == 'pass',]; dim(openclinica)


# Load Mutect WES data -----------------------------------------------------

# read tsv file
df_duke <- readMutect(paste0(wes_datapath, "Duke_mutect_pindel/Duke/mutect"), source = 'Duke', pattern = openclinica$sample_name_pdcis[openclinica$batch == 'Duke1'])
df_nki <- readMutect(paste0(wes_datapath, "NKI_mutect_pindel/NKI/mutect"), source = 'NKI', pattern = openclinica$sample_name_pdcis[openclinica$batch %in% c('NKI1', 'NKI2', 'NKI3', 'NKI4')])
df_sloane <- readMutect(paste0(wes_datapath, "SLO_mutect_pindel/SLO/mutect"), source = 'SLO', pattern = openclinica$sample_name_pdcis[openclinica$batch %in% c('SLO1', 'SLO2')])
df_melbourne <- readMutect(paste0(wes_datapath, "SLO3_mutect_pindel/SLO3/mutect"), source = 'SLO3', pattern = openclinica$sample_name_pdcis[openclinica$batch == 'SLO3'])

# merge data
df_mrg <- rbind(df_duke, df_nki, df_sloane, df_melbourne); dim(df_mrg)
ifelse(length(unique(openclinica$sampleid_pdcis)) == length(unique(df_mrg$tumor_name)), 'Everything is okay', stop("There are samples which weren't read"))

# add clinical data
df_mrg <- merge(df_mrg, openclinica, by.x = 'tumor_name', by.y = 'sampleid_pdcis'); dim(df_mrg)
df_mrg$location <- paste0('chr', df_mrg$chr, ':', df_mrg$start, '-', df_mrg$end)


# Filter somatic mutations ------------------------------------------------

# add if hotspot 
dat_hotspots <- read.csv("./data/hotspots.csv") 
df_mrg$hotspot <- 'FALSE'
for (i in 1:nrow(dat_hotspots)){
  message(i)
  Ind <- which(df_mrg$gene.knowngene == dat_hotspots$gene[i] & df_mrg$aaannotation == dat_hotspots$AA_mut[i] & df_mrg$ref_allele == dat_hotspots$ref[i] & df_mrg$alt_allele == dat_hotspots$alt[i] & df_mrg$location == dat_hotspots$hg19[i])
  
  if(length(Ind) > 0){
    df_mrg$hotspot[Ind] <- 'TRUE'
  }
}

# annotate variants
#df_mrg$variant_classification = recode_variant_classification(df_mrg$exonicfunc.knowngene)
df_mrg <- df_mrg %>% mutate(key = paste0(chr, ":", start, ":", ref_allele, "-", alt_allele)) %>% 
  dplyr::mutate(aaannotation1 = gsub("[p[:punct:]]", "", aaannotation, ignore.case = FALSE),
                aaannotation2 = gsub("([A-Z]?)([0-9]{,6})(.*)", "\\1\\2", aaannotation1),
                hotkey = paste0(gene.knowngene, "-", aaannotation2),
                normal_af = n_alt_count/(n_ref_count + n_alt_count),
                t_depth = t_ref_count + t_alt_count,
                n_depth = n_ref_count + n_alt_count)

is.na(df_mrg$key) %>% table()

# filter variants in gene panel
df_mrg <- filterByBed(anned=df_mrg, bed=bed, chrom = 'chr', start = 'start', end = 'end')

# remove mutations with low coverage
df_mrg <- df_mrg[-which(MyNumeric(df_mrg$t_depth) < 15 & df_mrg$hotspot == 'FALSE'),]
df_mrg <- df_mrg[-which(MyNumeric(df_mrg$n_depth) < 10 & df_mrg$hotspot == 'FALSE'),]
df_mrg <- df_mrg[-which(MyNumeric(df_mrg$normal_af) >= 0.01 & df_mrg$hotspot == 'FALSE'),]

# removing synonymous snps
df_mrg <- df_mrg[!df_mrg$exonicfunc.knowngene %in% c("synonymous SNV", 'unknown'),]

# filtering exonic and splicing
df_mrg <- df_mrg[df_mrg$func.knowngene %in% c('exonic', 'splicing', 'exonic;splicing'),]
df_mrg$func.knowngene <- ifelse(df_mrg$func.knowngene == 'exonic;splicing', 'exonic', df_mrg$func.knowngene)  # "splicing" in ANNOVAR is defined as variant that is within 2-bp away from an exon/intron boundary by default, but the threshold can be changed by the --splicing_threshold argument. Before Feb 2013, if "exonic,splicing" is shown, it means that this is a variant within exon but close to exon/intron boundary; this behavior is due to historical reason, when a user requested that exonic variants near splicing sites be annotated with splicing as well. However, I continue to get user emails complaining about this behavior despite my best efforts to put explanation in the ANNOVAR website with details. Therefore, starting from Feb 2013 , "splicing" only refers to the 2bp in the intron that is close to an exon, and if you want to have the same behavior as before, add -exonicsplicing argument.

# remove mutations with less than 10 log odds score
df_mrg <- df_mrg[-which(MyNumeric(df_mrg$t_lod_fstar) < 10 & df_mrg$hotspot == 'FALSE'),]

# remove mutations in esp6500 database
df_mrg <- df_mrg[-which(MyNumeric(df_mrg$esp6500siv2_all) >= 0.01 & df_mrg$hotspot == 'FALSE'),]

# remove mutations in exac database
df_mrg <- df_mrg[-which(MyNumeric(df_mrg$exac_all) >= 0.01 & df_mrg$hotspot == 'FALSE'),]

# remove mutations in 100G database
df_mrg <- df_mrg[-which(MyNumeric(df_mrg$x1kg2015aug_max) >= 0.01 & df_mrg$hotspot == 'FALSE'),]

# remove mutations with less than 2% vafs 
df_mrg <- df_mrg[-which(as.numeric(df_mrg$tumor_f) < 0.02 & df_mrg$hotspot == 'FALSE'),]; dim(df_mrg)

# filter mutations in cosmic
df_mrg <- df_mrg[!(df_mrg$ref_allele=="C" & df_mrg$alt_allele=="T" & GetCosmicNumber(df_mrg$cosmic70)<3 & as.numeric(df_mrg$tumor_f)<0.1 & df_mrg$hotspot == 'FALSE'),]; dim(df_mrg)
df_mrg <- df_mrg[!(df_mrg$ref_allele=="G" & df_mrg$alt_allele=="A" & GetCosmicNumber(df_mrg$cosmic70)<3 & as.numeric(df_mrg$tumor_f)<0.1 & df_mrg$hotspot == 'FALSE'),]; dim(df_mrg)

# keep splice variants if +5 of gene start and end
grch37 = useEnsembl(biomart="genes",GRCh=37, dataset="hsapiens_gene_ensembl")
gene_annotation <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                        'start_position', 'end_position', 'band','external_synonym'),
                         filters = 'hgnc_symbol', 
                         values = unique(bed[,4]), 
                         mart = grch37)
gene_annotation_2 <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                          'start_position', 'end_position', 'band','external_synonym'),
                           filters = 'external_synonym', 
                           values = unique(bed$V4), 
                           mart = grch37)
gene_annotation <- rbind(gene_annotation, gene_annotation_2)
gene_annotation <- rbind(gene_annotation, data.frame(hgnc_symbol = 'MRTFA', chromosome_name = 'chr22', start_position=40806285, end_position=41032723, band=NA, external_synonym='MKL1'))
gene_annotation$hgnc_symbol <- ifelse(gene_annotation$hgnc_symbol=='MLLT4', 'AFDN', gene_annotation$hgnc_symbol)
gene_annotation <- gene_annotation[gene_annotation$chromosome_name %in% c(1:22, 'X'),]

Ind <- which(df_mrg$func.knowngene == 'splicing')
to_remove <- c()
for(i in Ind){
  gene <- strsplit(df_mrg$gene[i], ';')[[1]][1]
  query <- gene_annotation[gene_annotation$hgnc_symbol==gene,]
  if(nrow(query)==0){
    query <- gene_annotation[gene_annotation$external_synonym==gene,]
  } else if (nrow(query) > 1){
    query <- query[1,]
  }
  if (!isTRUE(df_mrg$start[i] > query$start_position+5 & df_mrg$end[i] < query$end_position-5)){
    to_remove <- c(to_remove,i)
  }
}

if(!is.null(to_remove)){
  df <- df[-to_remove,]
}


# export filtered mutations
saveRDS(df_mrg, "./data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds")
write.table(df_mrg, './data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# Load Pindel WES data -----------------------------------------------------

# read tsv file
df_duke <- readMutect(paste0(wes_datapath, "Duke_mutect_pindel/Duke/pindel"), source = 'Duke', pattern = openclinica$sample_name_pdcis[openclinica$batch == 'Duke1'])
df_nki <- readMutect(paste0(wes_datapath, "NKI_mutect_pindel/NKI/pindel"), source = 'NKI', pattern = openclinica$sample_name_pdcis[openclinica$batch %in% c('NKI1', 'NKI2', 'NKI3', 'NKI4')])
df_sloane <- readMutect(paste0(wes_datapath, "SLO_mutect_pindel/SLO/pindel"), source = 'SLO', pattern = openclinica$sample_name_pdcis[openclinica$batch %in% c('SLO1', 'SLO2')])
df_melbourne <- readMutect(paste0(wes_datapath, "SLO3_mutect_pindel/SLO3/pindel"), source = 'SLO3', pattern = openclinica$sample_name_pdcis[openclinica$batch == 'SLO3'])

# merge data
df_indels <- rbind(df_duke, df_nki, df_sloane, df_melbourne); dim(df_indels)
df_indels$tumor_name <- gsub('Futreal-Precision-clonalrela-','',gsub('AFutreal-PRECISION-Heterogeneity-','',gsub('ESawyer-PrecDCIS-DNA-','',gsub('Lips-PRECISION-NKIWES2-', '', gsub('Lips-PRECISION-NKIWES3-', '', gsub('Lips-PRECISION-NKIWES2-','',gsub('LEsther-PRECISION-NKIWES4-', '', gsub('DShelleyHwang-PRECISION-WES-' ,'' , df_indels$sample_name))))))))
ifelse(length(unique(openclinica$sampleid_pdcis)) == length(unique(df_indels$tumor_name)), 'Everything is okay', stop("There are samples which weren't read"))

# add clinical data
df_indels2 <- merge(df_indels, openclinica, by.x = 'tumor_name', by.y = 'sampleid_pdcis'); dim(df_indels2)

# annotate variants
df_indels2 <- df_indels2 %>% mutate(key = paste0(chrom, ":", start, "-", end, ":", type, "-", length)) %>% dplyr::mutate(aaannotation1 = gsub("[p[:punct:]]", "", aaannotation, ignore.case = FALSE),
                                                                                                                       aaannotation2 = gsub("([A-Z]?)([0-9]{,6})(.*)", "\\1\\2", aaannotation1),
                                                                                                                       hotkey = paste0(gene.knowngene, "-", aaannotation2),
                                                                                                                       normal_af = n_alt_count/(n_ref_count + n_alt_count),
                                                                                                                       t_depth = t_ref_count + t_alt_count,
                                                                                                                       n_depth = n_ref_count + n_alt_count)

df_indels2$matching_key = paste0(df_indels2$tumor_name, "_", df_indels2$key, "_", df_indels2$hotkey)

# filtering
df_indels2 <- df_indels2 %>% dplyr::filter(MyNumeric(t_vaf) >= 0.02, #remove allele fraction in tumor samples < 0.02
                                           MyNumeric(esp6500siv2_all) < 0.01, # remove mutations in esp6500 database
                                           MyNumeric(exac_all) < 0.01,  # remove mutations in exac databases
                                           MyNumeric(x1kg2015aug_max) < 0.01); dim(df_indels2) # remove mutations in 100G database

# remove mutations with low coverage
df_indels2 <- df_indels2 %>% dplyr::filter(MyNumeric(t_depth) >= 15, MyNumeric(n_depth) >= 10, MyNumeric(n_vaf) <= 0.01); dim(df_indels2) 

# keep exonic mutations
df_indels2 <- df_indels2 %>% dplyr::filter(func.knowngene %in% c("exonic")); dim(df_indels2) 

# remove tandem duplications
df_indels_notds <- df_indels2 %>% dplyr::filter(!type %in% c("TD")); dim(df_indels2) 
df_indels_notds = df_indels_notds %>% dplyr::group_by(patient_id) %>% mutate(mut_cnt = n())

# filter for strand bias
df_indels_notds1 = df_indels_notds %>% dplyr::filter(supportup >= 2, supportdown >= 2) 

# filter non synonymous mutations
df_indels_notds1 <- df_indels_notds1[!df_indels_notds1$exonicfunc.knowngene %in% c("synonymous SNV", 'unknown'),]

# export filtered mutations
saveRDS(df_indels_notds1, "./data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds")
write.table(df_indels_notds1, './data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)




