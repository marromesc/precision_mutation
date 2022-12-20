# load libraries ----------------------------------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, reshape2)
p_load(glue, janitor, params)
p_load(limma, broom)
p_load(forcats, effsize)
p_load(MultiAssayExperiment, SummarizedExperiment)
p_load(ComplexHeatmap, cowplot, ggsci, ggpubr, ggstatsplot)
p_load(LBSPR, magrittr,maftools, stargazer)
library(gridExtra);library(grid);library(tidylog)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
source("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/scripts/00_functions.R")
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/")
dir.create("01_read_files_v2")
odir = "01_read_files_v2"
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/01_read_files_v2/")

source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/mut2maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/to_vcf.maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_sheets.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_mutect.R')
#install.packages('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/github_wranglr/', repos = NULL, type="source")
library(wranglr)


# path to files ----------------------------------------------------------------------------------------------------
path_to_duke = "~/sea_rsrch1/sandbox/20200622_dcis_exome/mutect/"
path_to_duke_sloane = "~/sea_rsrch1/sandbox/20200216_dcis_trios_batch2/mutect/"
path_to_sloane = "~/sea_rsrch1/sandbox/dcis_exome/kings/mutect/"
path_to_nki = "~/sea_rsrch1/sandbox/20200110_dcis_trios/mutect/"
path_to_nki2 = "~/sea_rsrch1/sandbox/20190520_dcis_exome/190515_A00482_0065_AHJ3YVDMXX/mutect/"
path_to_mda = "~/sea_rsrch1/data/dcis/190110_A00482_0045_AHGJCGDMXX/mutect"

path_to_nki3 = "~/sea_rsrch1/sandbox/210909_Futreal_LipsEsther_PRECISION_WES_box/mutect/"
path_to_nki4 = "~/sea_rsrch1/sandbox/211116_Futreal_LipsEsther_PRECISION_WES_box/mutect/"
path_to_sloane2 = "~/sea_rsrch1/sandbox/210909_ESawyer_PrecisionDCIS_WES_box/mutect/"
path_to_sloane3 = "~/sea_rsrch1/sandbox/211116_ESawyer_PrecisionDCIS_WES_box/mutect/"

# read files ----------------------------------------------------------------------------------------------------
df_mda = read_mutect(datapath = path_to_mda, source = "mda")
df_duke = read_mutect(datapath = path_to_duke, source = "duke1")
df_duke_sloane = read_mutect(datapath = path_to_duke_sloane, source = "duke_sloane")
df_sloane = read_mutect(datapath = path_to_sloane, source = "sloane")
df_nki = read_mutect(datapath = path_to_nki, source = "nki1")
df_nki2 = read_mutect(datapath = path_to_nki2, source = "nki2")

df_nki3 = read_mutect(datapath = path_to_nki3, source = "nki3")
df_nki4 = read_mutect(datapath = path_to_nki4, source = "nki4")
df_sloane2 = read_mutect(datapath = path_to_sloane2, source = "sloane2")
df_sloane3 = read_mutect(datapath = path_to_sloane3, source = "sloane3")

# annotate sample names
df_mda = df_mda %>% mutate(sample_idx = str_extract(sample_name, "PRE3.*"))
df_duke = df_duke %>% mutate(sample_idx = str_extract(sample_name, "DCIS.*"))
df_nki2 = df_nki2 %>% mutate(sample_idx = str_extract(sample_name, "NKI.*"))
df_duke_sloane = df_duke_sloane %>% mutate(sample_idx = str_extract(sample_name, "WES-.*|S111.*"))
df_sloane = df_sloane %>% mutate(sample_idx = str_extract(sample_name, "clonalrela-.*"))

df_nki = df_nki %>% mutate(sample_idx = str_extract(sample_name, "PRE-.*"))
df_nki = df_nki %>% mutate(sample_idx = str_extract(sample_idx, "NKI.*"))

df_nki3 = df_nki3 %>% mutate(sample_idx = gsub("LEsther-PRECISION-NKIWES4-PRE-NKI-(.*)", "\\1", sample_name))
df_nki4 = df_nki4 %>% mutate(sample_idx = gsub("LEsther-PRECISION-NKIWES4-PRE-NKI-(.*)", "\\1", sample_name))

df_sloane2 = df_sloane2 %>% mutate(sample_idx = gsub("ESawyer-PrecDCIS-DNA-(.*)", "\\1", sample_name))
df_sloane3 = df_sloane3 %>% mutate(sample_idx = gsub("ESawyer-PrecDCIS-DNA-(.*)", "\\1", sample_name))

# gsub give pattern i.e common, then (.*) matches anything., \\1 extract 1st part
# f = df_nki3$sample_name[1]
# gsub("LEsther-PRECISION-NKIWES4-PRE-NKI-(.*)", "\\1", f)


# merge data
df_mrg = bind_rows(df_mda, df_duke)
df_mrg = bind_rows(df_mrg, df_nki2)
df_mrg = bind_rows(df_mrg, df_nki)
df_mrg = bind_rows(df_mrg, df_duke_sloane)
df_mrg = bind_rows(df_mrg, df_sloane)

df_mrg = bind_rows(df_mrg, df_sloane2)
df_mrg = bind_rows(df_mrg, df_sloane3)
df_mrg = bind_rows(df_mrg, df_nki3)
df_mrg = bind_rows(df_mrg, df_nki4)
write_rds(df_mrg, "df_mrg.rds") #358920    162, #378 samples
write.csv(df_mrg, "df_mrg.csv")
write.csv(data.frame(sample_idx = unique(df_mrg$sample_idx)), "sampleidx.csv")
# add source manually (dumb coding)
sample_idx = read.csv("sampleidx.csv")
sample_idx = sample_idx[, -1]
head(sample_idx)

# filter data -------------------------------------------------------------------------------------------
df_mrg1 = left_join(df_mrg, sample_idx, by = "sample_idx") #358920    167
df_mrg2 = df_mrg1 %>% dplyr::filter(include == "yes", 
                                    case_cntrl %in% c("case", "control"), 
                                    timepoint == "primary")
dim(df_mrg2) #206780    167
write_rds(df_mrg2, "df_mrg2.rds")
write_rds(df_mrg2, "df_mrg2.csv")

# keep sample that pass QC
df_mrg2 =  df_mrg2 %>% dplyr::filter(!sample_idx %in% c("PRE30006-03-sd-01-T", "PRE30013-03-sd-01-T",
                                                         "PRE30004-03-sd-01-T", "PRE30009-03-sd-01-T",
                                                         "PRE30012-03-sd-01-T", "PRE30011-03-sd-01-T",
                                                         "NKI92-PRI-T", "NKI9117-T", "NKI8341-PRI-T",
                                                         "NKI8143-PRI-T", "NKI809-PRI-T",
                                                         "NKI7831-T", "NKI783-T", "NKI7788-PRI-T",
                                                         "NKI7472-PRI-T", "NKI63-T", "NKI5753-PRI-T",
                                                         "NKI5320-PRI-T", "NKI4410-PRI-T", "NKI3970-PRI-T",
                                                         "NKI2406-PRI-T", "NKI2354-PRI-T",
                                                         "NKI10067-PRI-T", "S111113-T", "S111135-T",
                                                         "S111585-T", "S111531-T", "NKI7340-PRI-T",
                                                         "NKI3669-PRI-T", "NKI3409-PRI-T", "NKI2175-PRI-T",
                                                         "NKI1664-PRI-T", "NKI1292-PRI-T", "NKI10126-PRI-T",
                                                         "DCIS00005Primary-T", "DCIS00511Primary-T", "DCIS00513Primary-T",
                                                         "DCIS00509Primary-T", "DCIS00501Primary-T", "DCIS00502Primary-T",
                                                         "DCIS00496Primary-T", "DCIS00493Primary-T", "DCIS00492Primary-T",
                                                         "DCIS00481Primary-T", "DCIS00402Primary-T", "DCIS00309Primary-T",
                                                         "DCIS00305Primary-T", "DCIS00276Primary-T", "DCIS00219Primary-T", 
                                                         "DCIS00162Primary-T", "DCIS00142Primary-T", "DCIS00135Primary-T",
                                                         "DCIS00047Primary-T", "DCIS00005Primary-T"))

write_rds(df_mrg2, "df_mrg2.rds")
dim(df_mrg2) #151915    167
length(unique(df_mrg2$sample_idx)) #176
# 119 cases, 57 controls
# 29 DCIS cases, 90 inv cases, 

is.na(df_mrg2$sample_idx) %>% table()
is.na(df_mrg2$case_cntrl) %>% table()
is.na(df_mrg2$recurrence) %>% table()


# annotate variants
df_mrg2$variant_classification = recode_variant_classification(df_mrg2$exonicfunc.knowngene)
df_mrg2 = df_mrg2 %>% mutate(key = paste0(chr, ":", start, ":", ref_allele, "-", alt_allele)) %>% 
        dplyr::mutate(aaannotation1 = gsub("[p[:punct:]]", "", aaannotation, ignore.case = FALSE),
                      aaannotation2 = gsub("([A-Z]?)([0-9]{,6})(.*)", "\\1\\2", aaannotation1),
                      hotkey = paste0(gene.knowngene, "-", aaannotation2),
                      normal_af = n_alt_count/(n_ref_count + n_alt_count),
                      t_depth = t_ref_count + t_alt_count,
                      n_depth = n_ref_count + n_alt_count)

is.na(df_mrg2$key) %>% table()

# filter data -------------------------------------------------------------------------------------------
df_mrg3 = df_mrg2 %>% dplyr::filter(exonicfunc.knowngene %in% c("synonymous SNV","nonsynonymous SNV", "stoploss", "stopgain"))
dim(df_mrg3) #118059 172

df_mrg3$exac_all = as.numeric(df_mrg3$exac_all)
df_mrg3$t_lod_fstar = as.numeric(df_mrg3$t_lod_fstar)
df_mrg3$esp6500siv2_all = as.numeric(df_mrg3$esp6500siv2_all)
df_mrg3$x1kg2015aug_max = as.numeric(df_mrg3$x1kg2015aug_max)
df_mrg3$t_depth = as.numeric(df_mrg3$t_depth)
df_mrg3$n_depth = as.numeric(df_mrg3$n_depth)
df_mrg3$normal_af = as.numeric(df_mrg3$normal_af)
df_mrg3$tumor_f = as.numeric(df_mrg3$tumor_f)

# remove samples with less than 2% vafs 
df_mrg4 = df_mrg3 %>% dplyr::filter(t_lod_fstar >= 10);dim(df_mrg3);dim(df_mrg4)

# remove samples in exac databases 
df_mrg5 = df_mrg4 %>% dplyr::filter(tumor_f >= 0.02, esp6500siv2_all < 0.01 | is.na(esp6500siv2_all), 
                                    exac_all < 0.01 | is.na(exac_all), 
                                    x1kg2015aug_max < 0.01 | is.na(x1kg2015aug_max))
df_mrg6 = df_mrg5 %>% dplyr::filter(t_depth >= 15, n_depth >= 10, normal_af < 0.01)

dim(df_mrg5)

write_rds(df_mrg5, "df_mrg5.rds")


