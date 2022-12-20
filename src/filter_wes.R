library(tidyr)
library(dplyr)
library(stringr)
library(gridExtra)

read_mutect = function(datapath, source, pattern = NULL){
  # combine all mutects 
  if (is.null(pattern)){
    ls.files = list.files(path = c(datapath), pattern = '.tsv', full.names = TRUE)
  } else {
    ls.files = do.call(rbind, lapply(pattern, function(x) list.files(path = c(datapath),
                                                                     pattern = x, full.names = TRUE)))[,1]
  }
  
  datlist = list()
  for(i in 1:length(ls.files)){
    #message(ls.files[i])
    if( file.info(ls.files[i])$size > 0 ){
      dat.mat = read.delim2(file = ls.files[i]) %>% 
        tbl_df() %>% 
        mutate(tool = "mutect", source = source)
      
      datlist[[i]] = dat.mat
    } 
  }
  
  df_comb = do.call(rbind, datlist);dim(df_comb) 
  return(df_comb)
}

setwd('/mnt/albyn/maria/precision_mutation')

# read openclinica master sheet
openclinica <- as.data.frame(readxl::read_xlsx('./data/updated_20221114_clinicaldata_caco_genomic_samples.xlsx')); dim(openclinica)
openclinica2 <- as.data.frame(readxl::read_xlsx('./data/updated_20221114_clinicaldata_caco_genomic_samples.xlsx', sheet = 2))
openclinica2 <- openclinica2[openclinica2$redcap_event == 'PRI',]
openclinica <- as.data.frame(merge(openclinica, openclinica2)); dim(openclinica)


# Perform QC --------------------------------------------------------------

# read qc master sheet
sloane_qc <- readxl::read_xlsx('./data/WES/20221102_DCIS_Precision_WES_sampleInfo_by_XiaogangWu.xlsx', sheet = 'SLO')
nki_qc <- readxl::read_xlsx('./data/WES/20221102_DCIS_Precision_WES_sampleInfo_by_XiaogangWu.xlsx', sheet = 'NKI')
nki_qc <- nki_qc[!(nki_qc$note %in% c("No sample of primary DCIS, only DCIS recurrence.", "\"not usable\"")),] # remove not usable samples
duke_qc <- readxl::read_xlsx('./data/WES/20221102_DCIS_Precision_WES_sampleInfo_by_XiaogangWu.xlsx', sheet = 'Duke')

wes_qc <- rbind(sloane_qc[,1:20], nki_qc[,1:20], duke_qc[,1:20])

# Perform QC (Total reads >= 10M AND median target coverage >= 40 for normal and >= 80 for tumor samples)
wes_qc$qc <- ifelse(wes_qc$`total_reads(million)` >= 10 & wes_qc$tissue_pathology == 'Normal' & wes_qc$median_Target_Cov >= 40, 'pass', 
                ifelse(wes_qc$`total_reads(million)` >= 10 & wes_qc$tissue_pathology != 'Normal' & wes_qc$median_Target_Cov >= 80, 'pass', 'fail'))

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
openclinica <- openclinica[openclinica$first_subseq_event %in% c('ipsilateral DCIS', 'ipsilateral IBC', 'NA', 'death'),]; dim(openclinica)

# filter samples which passed qc
openclinica <- openclinica[openclinica$qc_normal == 'pass' & openclinica$qc_pdcis == 'pass',]; dim(openclinica)


# Load Mutect WES data -----------------------------------------------------

# read tsv file
df_duke <- read_mutect("./data/WES/DCIS_Precision_WES_Duke/mutect", source = 'Duke', pattern = openclinica$sample_name_pdcis[openclinica$batch == 'Duke1'])
df_nki <- read_mutect("./data/WES/DCIS_Precision_WES_NKI/mutect", source = 'NKI', pattern = openclinica$sample_name_pdcis[openclinica$batch %in% c('NKI1', 'NKI2', 'NKI3', 'NKI4')])
df_sloane <- read_mutect("./data/WES/DCIS_Precision_WES_SLO/mutect", source = 'SLO', pattern = openclinica$sample_name_pdcis[openclinica$batch %in% c('SLO1', 'SLO2')])
df_melbourne <- read_mutect("./data/WES/DCIS_Precision_WES_SLO3/mutect", source = 'SLO3', pattern = openclinica$sample_name_pdcis[openclinica$batch == 'SLO3'])

# merge data
df_mrg <- rbind(df_duke, df_nki, df_sloane, df_melbourne); dim(df_mrg)
ifelse(length(unique(openclinica$sampleid_pdcis)) == length(unique(df_mrg$tumor_name)), 'Everything is okay', stop("There are samples which weren't read"))

# add clinical data
df_mrg <- merge(df_mrg, openclinica, by.x = 'tumor_name', by.y = 'sampleid_pdcis'); dim(df_mrg)


# Filter somatic mutations ------------------------------------------------

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

df_mrg$exac_all <- as.numeric(df_mrg$exac_all)
df_mrg$t_lod_fstar <- as.numeric(df_mrg$t_lod_fstar)
df_mrg$esp6500siv2_all <- as.numeric(df_mrg$esp6500siv2_all)
df_mrg$x1kg2015aug_max <- as.numeric(df_mrg$x1kg2015aug_max)
df_mrg$t_depth <- as.numeric(df_mrg$t_depth)
df_mrg$n_depth <- as.numeric(df_mrg$n_depth)
df_mrg$normal_af <- as.numeric(df_mrg$normal_af)
df_mrg$tumor_f <- as.numeric(df_mrg$tumor_f)

# save merged mutations
write.table(df_mrg, './data/WES/DCIS_Precision_CaCo_Panel_Sloane_Mutect.txt', sep = '\t', quote = FALSE, row.names = FALSE)

# keep snp mutations 
df_mrg2 <- df_mrg %>% dplyr::filter(func.knowngene %in% c('exonic', 'splicing', 'exonic;splicing')); dim(df_mrg2); unique(df_mrg2$exonicfunc.knowngene)

# keep snp mutations 
#df_mrg3 <- df_mrg2 %>% dplyr::filter(exonicfunc.knowngene %in% c('.', "synonymous SNV","nonsynonymous SNV", "stoploss", "stopgain")); dim(df_mrg3); unique(df_mrg3$exonicfunc.knowngene)

# remove mutations with less than 10 log odds score
df_mrg4 <- df_mrg2 %>% dplyr::filter(t_lod_fstar >= 10); dim(df_mrg4) ; unique(df_mrg4$exonicfunc.knowngene)

df_mrg5 <- df_mrg4 %>% dplyr::filter(tumor_f >= 0.02, # remove mutations with less than 2% vafs 
                                     esp6500siv2_all < 0.01 | is.na(esp6500siv2_all),  # remove mutations in esp6500 database
                                     exac_all < 0.01 | is.na(exac_all),  # remove mutations in exac database
                                     x1kg2015aug_max < 0.01 | is.na(x1kg2015aug_max)); dim(df_mrg5); unique(df_mrg5$exonicfunc.knowngene) # remove mutations in 100G database
df_mrg6 <- df_mrg5 %>% dplyr::filter(t_depth >= 15, n_depth >= 10, normal_af < 0.01); dim(df_mrg6); unique(df_mrg6$exonicfunc.knowngene) # remove mutations with low coverage

# filter non synonymous mutations
df_mrg6 <- df_mrg6[!df_mrg6$exonicfunc.knowngene %in% c("synonymous SNV", 'unknown'),]; dim(df_mrg6); unique(df_mrg6$exonicfunc.knowngene)

# rename mutations
df_mrg6$exonicfunc.knowngene[df_mrg6$exonicfunc.knowngene %in% c(".")] <- "splicing"
df_mrg6$exonicfunc.knowngene[df_mrg6$exonicfunc.knowngene %in% c("nonsynonymous SNV")] <- "missense"
df_mrg6$exonicfunc.knowngene[df_mrg6$exonicfunc.knowngene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"
unique(df_mrg6$exonicfunc.knowngene)

# export filtered mutations
saveRDS(df_mrg6, "./data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds")
write.table(df_mrg6, './data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# Load Pindel WES data -----------------------------------------------------

# read tsv file
df_duke <- read_mutect("./data/WES/DCIS_Precision_WES_Duke/pindel", source = 'Duke', pattern = openclinica$sample_name_pdcis[openclinica$batch == 'Duke1'])
df_nki <- read_mutect("./data/WES/DCIS_Precision_WES_NKI/pindel", source = 'NKI', pattern = openclinica$sample_name_pdcis[openclinica$batch %in% c('NKI1', 'NKI2', 'NKI3', 'NKI4')])
df_sloane <- read_mutect("./data/WES/DCIS_Precision_WES_SLO/pindel", source = 'SLO', pattern = openclinica$sample_name_pdcis[openclinica$batch %in% c('SLO1', 'SLO2')])
df_melbourne <- read_mutect("./data/WES/DCIS_Precision_WES_SLO3/pindel", source = 'SLO3', pattern = openclinica$sample_name_pdcis[openclinica$batch == 'SLO3'])

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

df_indels2$exac_all <- as.numeric(df_indels2$exac_all)
df_indels2$esp6500siv2_all <- as.numeric(df_indels2$esp6500siv2_all)
df_indels2$x1kg2015aug_max <- as.numeric(df_indels2$x1kg2015aug_max)
df_indels2$t_depth <- as.numeric(df_indels2$t_depth)
df_indels2$n_depth <- as.numeric(df_indels2$n_depth)
df_indels2$normal_af <- as.numeric(df_indels2$normal_af)
df_indels2$tumor_f <- as.numeric(df_indels2$tumor_f)
df_indels2$t_vaf <- as.numeric(df_indels2$t_vaf)
df_indels2$n_vaf <- as.numeric(df_indels2$n_vaf)

df_indels2$matching_key = paste0(df_indels2$tumor_name, "_", df_indels2$key, "_", df_indels2$hotkey)

# filtering
df_indels2 <- df_indels2 %>% dplyr::filter(t_vaf >= 0.02, #remove allele fraction in tumor samples < 0.02
                                           esp6500siv2_all < 0.01 | is.na(esp6500siv2_all), # remove mutations in esp6500 database
                                           exac_all < 0.01 | is.na(exac_all),  # remove mutations in exac databases
                                           x1kg2015aug_max < 0.01 | is.na(x1kg2015aug_max)); dim(df_indels2) # remove mutations in 100G database

# remove mutations with low coverage
df_indels2 <- df_indels2 %>% dplyr::filter(t_depth >= 15, n_depth >= 10, n_vaf <= 0.01); dim(df_indels2) 

# keep exonic mutations
df_indels2 <- df_indels2 %>% dplyr::filter(func.knowngene %in% c("exonic")); dim(df_indels2) 

# remove tandem duplications
df_indels_notds <- df_indels2 %>% dplyr::filter(!type %in% c("TD")); dim(df_indels2) 
df_indels_notds = df_indels_notds %>% dplyr::group_by(patient_id) %>% mutate(mut_cnt = n())

# filter for strand bias
df_indels_notds1 = df_indels_notds %>% dplyr::filter(supportup >= 2, supportdown >= 2) 

# filter non synonymous mutations
df_indels_notds1 <- df_indels_notds1[!df_indels_notds1$exonicfunc.knowngene %in% c("synonymous SNV", 'unknown'),]

# rename mutations
df_indels_notds1$exonicfunc.knowngene[df_indels_notds1$exonicfunc.knowngene %in% c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution")] <- "inframe_indel"
df_indels_notds1$exonicfunc.knowngene[df_indels_notds1$exonicfunc.knowngene %in% c("frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift"
df_indels_notds1$exonicfunc.knowngene[df_indels_notds1$exonicfunc.knowngene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"

# export filtered mutations
saveRDS(df_indels_notds1, "./data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds")
write.table(df_indels_notds1, './data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# hotspots 
# dat_hotspots = read_delim(file = "~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/data/hotspots.txt", delim = "\t")  %>% as_tibble();dim(dat_hotspots) # 1165    7
# dat_hotspots_breast = dat_hotspots %>% clean_names() %>% 
#   dplyr::filter(str_detect(tumor_type_composition, "breast"), type == "single residue") %>% 
#   dplyr::mutate(hotkey = paste0(gene, "-", residue))
# write.csv(dat_hotspots_breast, "dat_hotspots_breast.csv")
# 
# hotspots = unique(dat_hotspots_breast$hotkey);length(hotspots) #385
# head(dat_hotspots_breast)


