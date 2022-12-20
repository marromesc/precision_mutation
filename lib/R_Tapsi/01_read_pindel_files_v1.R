# load libraries ----------------------------------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, reshape2, dplyr)
p_load(glue, janitor, params)
p_load(limma, broom)
p_load(forcats, effsize)
p_load(MultiAssayExperiment, SummarizedExperiment)
p_load(ComplexHeatmap, cowplot, ggsci, ggpubr, ggstatsplot)
p_load(LBSPR, magrittr,maftools, stargazer)
library(gridExtra);library(grid)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/analysis_indels/")
dir.create("01_read_pindel_files_v1")
odir = "01_read_pindel_files_v1"
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/analysis_indels/01_read_pindel_files_v1/")

source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/mut2maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/to_vcf.maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_sheets.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_mutect.R')
install.packages('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/github_wranglr/', repos = NULL, type="source")
library(wranglr)


# sample_tracker ---------
df_tracker = read.csv("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/analysis_batch1/sample_tracker.csv", sep = "\t")


# hotspots -------------
dat_hotspots = read_delim(file = "~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/data/hotspots.txt", delim = "\t")  %>% as_tibble();dim(dat_hotspots) # 1165    7
dat_hotspots_breast = dat_hotspots %>% clean_names() %>% 
        dplyr::filter(str_detect(tumor_type_composition, "breast"), type == "single residue") %>% 
        dplyr::mutate(hotkey = paste0(gene, "-", residue))
write.csv(dat_hotspots_breast, "dat_hotspots_breast.csv")

hotspots = unique(dat_hotspots_breast$hotkey);length(hotspots) #385
head(dat_hotspots_breast)

# get final samples -------
df_all_muts = read.csv("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/analysis_combined/01_read_files_v4/combined.csv")
df_all_muts = df_all_muts[,-c(1:3)]
final_pt_list = unique(df_all_muts$patient)
final_sample_list = unique(df_all_muts$sample_name)

# ____________---------
# **read files ------

## 1a. pindel indels merged file ------
trios_pt = final_pt_list

indel_filepath = "~/sea_rsrch1/sandbox/pindel/merged_filtered_pindel.tsv"
df_indels = read.delim2(file = indel_filepath) %>% 
        tibble() %>% 
        mutate(tool = "pindel", source = "all")

df_indels1 = df_indels        

# only get trios
df_indels1 = df_indels %>% dplyr::filter(sample_name %in% final_sample_list)

# extract sample id
df_indels2 = df_indels1 %>% mutate(sample_id = sub(".*Heterogeneity-|.*NKIWES2-|.*NKIWES3-|.*WES-|.*DNA-|.*clonalrela-", "", df_indels1$sample_name))
#df_indels = df_indels %>% mutate(sample_id = qdapRegex::ex_between(df_indels$sample_name, "Futreal-Precision-clonalrela-", "-T"))

df_indels3 = left_join(df_indels2, df_tracker, by = "sample_id");dim(df_indels3) #431034   117

df_indels3 = df_indels3 %>% mutate(timepoint = case_when(sample_id %in% "3040420IR-T" ~ "recurrence",
                                                             sample_id %in% "3040420P-T" ~ "primary",
                                                             sample_id %in% "8090730IR-T" ~ "recurrence",
                                                             sample_id %in% "8090730P-T" ~ "primary",
                                                             sample_id %in% "9100205IR-T" ~ "recurrence",
                                                             sample_id %in% "9100205P-T" ~ "primary",
                                                             sample_id %in% "9100923IR-T" ~ "recurrence",
                                                             sample_id %in% "9100923P-T" ~ "primary",
                                                             sample_id %in% "IT02006INVRec-T" ~ "recurrence",
                                                             sample_id %in% "IT02006P-T" ~ "primary", TRUE ~ as.character(timepoint)))

df_indels3 = df_indels3 %>% mutate(patient = case_when(sample_id %in% "3040420IR-T" ~ "3040420",
                                                           sample_id %in% "3040420P-T" ~ "3040420",
                                                           sample_id %in% "8090730IR-T" ~ "8090730",
                                                           sample_id %in% "8090730P-T" ~ "8090730",
                                                           sample_id %in% "9100205IR-T" ~ "9100205",
                                                           sample_id %in% "9100205P-T" ~ "9100205",
                                                           sample_id %in% "9100923IR-T" ~ "9100923",
                                                           sample_id %in% "9100923P-T" ~ "9100923",
                                                           sample_id %in% "IT02006INVRec-T" ~ "IT02006",
                                                           sample_id %in% "IT02006P-T" ~ "IT02006", TRUE ~ as.character(patient)))


df_indels3 = df_indels3 %>% mutate(recurrence_type = case_when(sample_id %in% "3040420IR-T" ~ "invasive",
                                                                   sample_id %in% "3040420P-T" ~ "invasive",
                                                                   sample_id %in% "8090730IR-T" ~ "invasive",
                                                                   sample_id %in% "8090730P-T" ~ "invasive",
                                                                   sample_id %in% "9100205IR-T" ~ "invasive",
                                                                   sample_id %in% "9100205P-T" ~ "invasive",
                                                                   sample_id %in% "9100923IR-T" ~ "invasive",
                                                                   sample_id %in% "9100923P-T" ~ "invasive",
                                                                   sample_id %in% "IT02006INVRec-T" ~ "invasive",
                                                                   sample_id %in% "IT02006P-T" ~ "invasive", TRUE ~ as.character(recurrence_type)))


df_indels3 = df_indels3 %>% mutate(sample_id = case_when(sample_id %in% "3040420IR-T" ~ "3040420IR",
                                                         sample_id %in% "3040420P-T" ~ "3040420P",
                                                         sample_id %in% "8090730IR-T" ~ "8090730IR",
                                                         sample_id %in% "8090730P-T" ~ "8090730P",
                                                         sample_id %in% "9100205IR-T" ~ "9100205IR",
                                                         sample_id %in% "9100205P-T" ~ "9100205P",
                                                         sample_id %in% "9100923IR-T" ~ "9100923IR",
                                                         sample_id %in% "9100923P-T" ~ "9100923P",
                                                         sample_id %in% "IT02006INVRec-T" ~ "IT02006INVRec",
                                                         sample_id %in% "IT02006P-T" ~ "IT02006P", TRUE ~ as.character(sample_id)))

meta_data = read_sheet("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/analysis_kings_duke_set2_run1/metdadata.xlsx", header = TRUE)

df_indels4 = df_indels3 %>% mutate(sample_id = case_when(grepl("DCIS-", sample_id) ~ paste0(sample_id, "_"), TRUE ~ sample_id))

df_indels5 = left_join(df_indels4, meta_data, by = c("sample_id" = "Pattern"))

df_indels5 = df_indels5 %>% mutate(timepoint = case_when(Type %in% "INV" ~ "recurrence",
                                                           Type %in% "DCIS" ~ "primary", TRUE ~ timepoint))
df_indels5 = df_indels5 %>% mutate(recurrence_type = case_when(Type %in% "INV" ~ "invasive",
                                                         Type %in% "DCIS" ~ "invasive", TRUE ~ recurrence_type))

df_indels5 = dplyr::select(df_indels5, -c(Individual, Sample, Type, sample, median_Target_Cov))
dim(df_indels5) #431034    117


# write combined mutect results file
write.csv(df_indels5, "df_indels5.csv");dim(df_indels5) #[1] 431034    117
write_rds(df_indels5, "df_indels5.rds");dim(df_indels5) #[1] 431034    117



