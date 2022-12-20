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
dir.create("01_start_with_qc")
odir = "01_start_with_qc"
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/01_start_with_qc/")

source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/mut2maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/to_vcf.maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_sheets.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_mutect.R')
#install.packages('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/github_wranglr/', repos = NULL, type="source")
library(wranglr)


# path to files ----------------------------------------------------------------------------------------------------
# path_to_duke = "~/sea_rsrch1/sandbox/20200622_dcis_exome/coverage/117.E87.bams.cov_summary.txt/"
# path_to_duke_sloane = "~/sea_rsrch1/sandbox/20200216_dcis_trios_batch2/mutect/"
# path_to_sloane = "~/sea_rsrch1/sandbox/dcis_exome/kings/mutect/"
# path_to_nki = "~/sea_rsrch1/sandbox/20200110_dcis_trios/mutect/"
# path_to_nki2 = "~/sea_rsrch1/sandbox/20190520_dcis_exome/190515_A00482_0065_AHJ3YVDMXX/mutect/"
# path_to_mda = "~/sea_rsrch1/data/dcis/190110_A00482_0045_AHGJCGDMXX/mutect"

path_to_nki3 = "~/sea_rsrch1/sandbox/210909_Futreal_LipsEsther_PRECISION_WES_box/210830_A00482_0250_AHLHKWDSX2_96.E1.PRECISION.bams.cov_summary.txt"
path_to_nki4 = "~/sea_rsrch1/sandbox/211116_Futreal_LipsEsther_PRECISION_WES_box/211111_A00482_0262_AH3KJWDSX3_9.E1.LEsther.bams.cov_summary.txt"
path_to_sloane2 = "~/sea_rsrch1/sandbox/210909_ESawyer_PrecisionDCIS_WES_box/210901_A00482_0251_AHKHMJDRXY_20.E1.bams.cov_summary.txt"
path_to_sloane3 = "~/sea_rsrch1/sandbox/211116_ESawyer_PrecisionDCIS_WES_box/211111_A00482_0262_AH3KJWDSX3_84.E1.ESawyer.bams.cov_summary.txt"

# read files ----------------------------------------------------------------------------------------------------
# df_mda = read_mutect(datapath = path_to_mda, source = "mda")
# df_duke = read_mutect(datapath = path_to_duke, source = "duke1")
# df_duke_sloane = read_mutect(datapath = path_to_duke_sloane, source = "duke_sloane")
# df_sloane = read_mutect(datapath = path_to_sloane, source = "sloane")
# df_nki = read_mutect(datapath = path_to_nki, source = "nki1")
# df_nki2 = read_mutect(datapath = path_to_nki2, source = "nki2")

df_nki3 = read.csv2(path_to_nki3, sep = "\t")
df_nki4 = read.csv2(path_to_nki4, sep = "\t")
df_sloane2 = read.csv2(path_to_sloane2, sep = "\t")
df_sloane3 = read.csv2(path_to_sloane3, sep = "\t")

# annotate sample names
# df_mda = df_mda %>% mutate(sample_idx = str_extract(sample_name, "PRE3.*"))
# df_duke = df_duke %>% mutate(sample_idx = str_extract(sample_name, "DCIS.*"))
# df_nki2 = df_nki2 %>% mutate(sample_idx = str_extract(sample_name, "NKI.*"))
# df_duke_sloane = df_duke_sloane %>% mutate(sample_idx = str_extract(sample_name, "WES-.*|S111.*"))
# df_sloane = df_sloane %>% mutate(sample_idx = str_extract(sample_name, "clonalrela-.*"))
# 
# df_nki = df_nki %>% mutate(sample_idx = str_extract(sample_name, "PRE-.*"))
# df_nki = df_nki %>% mutate(sample_idx = str_extract(sample_idx, "NKI.*"))

df_nki3 = df_nki3 %>% mutate(sample_idx = gsub("PRE-NKI-(.*)", "\\1", ExternalID))
df_nki4 = df_nki4 %>% mutate(sample_idx = gsub("PRE-NKI-(.*)", "\\1", ExternalID))

df_sloane2 = df_sloane2 %>% mutate(sample_idx = ExternalID)
df_sloane3 = df_sloane3 %>% mutate(sample_idx = ExternalID)

# sample_idx = read.csv("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/01_read_files_v2/sampleidx.csv")
# sample_idx = sample_idx[, -1]
# head(sample_idx)

# annotate sample ids ------
df_nki3 =  df_nki3 %>% dplyr::select(median_Target_Cov, sample_idx) %>% 
        tidyr::separate(col = sample_idx, sep = "-", into = c("newId","type")) %>% mutate(median_Target_Cov = as.numeric(median_Target_Cov))
df_nki4 =  df_nki4 %>% dplyr::select(median_Target_Cov, sample_idx) %>% 
        tidyr::separate(col = sample_idx, sep = "-", into = c("newId","type"))%>% mutate(median_Target_Cov = as.numeric(median_Target_Cov))
df_sloane2 =  df_sloane2 %>% dplyr::select(median_Target_Cov, sample_idx) %>% 
        tidyr::separate(col = sample_idx, sep="(?=N|P|R)", extra="merge", into = c("newId","type"))%>% mutate(median_Target_Cov = as.numeric(median_Target_Cov))
df_sloane3 =  df_sloane3 %>% dplyr::select(median_Target_Cov, sample_idx) %>% 
        tidyr::separate(col = sample_idx, sep="(?=N|P|R)", extra="merge", into = c("newId","type"))%>% mutate(median_Target_Cov = as.numeric(median_Target_Cov))
df_sloane3[is.na(df_sloane3)] <- "DCIS"

# make QC plots -------
df_nki = rbind(df_nki3, df_nki4)
df_nki$median_Target_Cov = round(df_nki$median_Target_Cov, digits = 0)
df_nki =  df_nki %>% mutate(new_Median = case_when(type == "NOR" & median_Target_Cov<40 ~ "<40x",
                                                    type %in% c("PRI", "RC1") & median_Target_Cov<100 ~ "<100x",
                                                    TRUE ~ "PASS"))
write.csv(df_nki, "df_nki.csv")

# g1 = ggplot(df_nki3, aes(x = type, y = newId, fill = cut(median_Target_Cov, c(0, 40, 100, Inf)))) +
#         geom_tile() + 
#         #scale_fill_continuous() +
#         geom_text(aes(label = median_Target_Cov), color = "black", size = 2) +
#         #coord_fixed() + 
# scale_fill_manual(name = "median_Target_Cov",
#                    values = c("(0,40]" = "orange",
#                               "(40,100]" = "yellow",
#                               "(100, Inf]" = "white"),
#                    labels = c("<= 40x", "40 < median_Target_Cov <= 100", "> 100")) + 
#         theme_cowplot()
g1 = ggplot(df_nki, aes(x = type, y = newId, fill = new_Median)) +
        geom_tile(color = "grey50",
                  lwd = 0.25,
                  linetype = 1) + 
        #scale_fill_continuous() +
        geom_text(aes(label = median_Target_Cov), color = "black", size = 3.5) +
        #coord_fixed() + 
        scale_fill_manual(name = "median_Target_Cov",
                          values = c("<40x" = "orange",
                                     "<100x" = "orange",
                                     "PASS" = "white"),
                          labels = c("<= 40x", "<100X", "pass")) + 
        theme_pubr();g1
ggsave(filename = "nki34_qc.png", plot = g1, width = 4, height = 8)

# sloane
df_sl = rbind(df_sloane2, df_sloane3)
df_sl$median_Target_Cov = round(df_sl$median_Target_Cov, digits = 0)
write.csv(df_sl, "df_sl.csv")
df_sl= df_sl %>% mutate(type = case_when(type == "Normal " ~ "Normal",
                                         type == "RINV " ~ "RINV",
                                         TRUE ~ type))
df_sl =  df_sl %>% mutate(new_Median = case_when(type == "Normal" & median_Target_Cov<40 ~ "<40x",
                                                   type %in% c("Primary", "DCIS", "RINV") & median_Target_Cov<100 ~ "<100x",
                                                   TRUE ~ "PASS"))

g1 = ggplot(df_sl, aes(x = type, y = newId, fill = new_Median)) +
        geom_tile(color = "grey50",
                  lwd = 0.25,
                  linetype = 1) + 
        #scale_fill_continuous() +
        geom_text(aes(label = median_Target_Cov), color = "black", size = 3.5) +
        #coord_fixed() + 
        scale_fill_manual(name = "median_Target_Cov",
                          values = c("<40x" = "orange",
                                     "<100x" = "orange",
                                     "PASS" = "white"),
                          labels = c("<= 40x", "<100X", "pass")) + 
        theme_pubr();g1
ggsave(filename = "sloane_qc.png", plot = g1, width = 6, height = 8)

# samples that failes for cc
nki_samples = df_nki %>% group_by(newId) %>% dplyr::filter(new_Median == "PASS", type %in% c("PRI", "NOR")) %>% mutate(nc = n()) %>% dplyr::filter(nc == 2) %>% pull(newId) %>% unique()
write.csv(nki_samples, "nki_samples_pass.csv") #39
sl_samples = df_sl %>% group_by(newId) %>% dplyr::filter(new_Median == "PASS", type %in% c("Primary", "Normal")) %>% mutate(nc = n()) %>% dplyr::filter(nc == 2) %>% pull(newId) %>% unique()
write.csv(sl_samples, "sl_samples_pass.csv") #4
