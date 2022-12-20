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
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/")
dir.create("02_annotate_mutations_v2")
odir = "02_annotate_mutations_v2"
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/02_annotate_mutations_v2/")

source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/mut2maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/to_vcf.maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_sheets.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_mutect.R')
#install.packages('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/github_wranglr/', repos = NULL, type="source")
library(wranglr)

# read files ----------------------------------------------------------------------------------------------------
sample_idx = read.csv("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/01_read_files_v2/sampleidx.csv")
sample_idx = sample_idx[, -1]

df_mrg5 = read_rds("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/01_read_files_v2/df_mrg5.rds")


df3_case = df_mrg5 %>% dplyr::filter(case_cntrl == "case", recurrence == "dcis")
df3_caseinv = df_mrg5 %>% dplyr::filter(case_cntrl == "case", recurrence == "invasive")
df3_control = df_mrg5 %>% dplyr::filter(case_cntrl == "control")
df_maf_f3.2_t1t0 = tsv2maf(df_mrg6,
                           gene = "gene", 
                           func = "variant_classification",
                           sample_name = "sample_idx",
                           ref_name = "normal_name")

df_maf_f3.2_t1t0_case = tsv2maf(df3_case,
                                gene = "gene", 
                                func = "variant_classification",
                                sample_name = "sample_idx",
                                ref_name = "normal_name")

df_maf_f3.2_t1t0_caseinv = tsv2maf(df3_caseinv,
                                gene = "gene", 
                                func = "variant_classification",
                                sample_name = "sample_idx",
                                ref_name = "normal_name")

df_maf_f3.2_t1t0_control = tsv2maf(df3_control,
                                   gene = "gene", 
                                   func = "variant_classification",
                                   sample_name = "sample_idx",
                                   ref_name = "normal_name")

df_maf_f3.2_t1t0_case$Tumor_Sample_Barcode %>% unique()
df_maf_f3.2_t1t0$Hugo_Symbol %>% unique()
df_maf_f3.2_t1t0$Tumor_Sample_Barcode %>% n_distinct()

# sample_idx1 = sample_idx
# sample_idx1 = sample_idx1[,-1]
colnames(sample_idx) = c("Tumor_Sample_Barcode", "source", "include", "recurrence", "case_cntrl", "timepoint")
maf_f3.2_t1t0 = read.maf(df_maf_f3.2_t1t0, isTCGA = F, clinicalData = sample_idx)
maf_f3.2_t1t0_case= read.maf(df_maf_f3.2_t1t0_case, isTCGA = F, clinicalData = sample_idx)
maf_f3.2_t1t0_caseinv= read.maf(df_maf_f3.2_t1t0_caseinv, isTCGA = F, clinicalData = sample_idx)
maf_f3.2_t1t0_control = read.maf(df_maf_f3.2_t1t0_control, isTCGA = F, clinicalData = sample_idx)

# maftools::write.mafSummary(maf_f3.2_t1t0, "df_maf_f3.2_t1t0.maf")
# maftools::write.mafSummary(df_maf_f3.2_t1t0_pr, "df_maf_f3.2_t1t0_pr.maf")



maf_f3.2_t1t0@data = maf_f3.2_t1t0@data %>% dplyr::filter(!Hugo_Symbol %in% c("TTN", "MUC6", "MUC2", "MUC17", "MUC4"))
maf_f3.2_t1t0_case@data = maf_f3.2_t1t0_case@data %>% dplyr::filter(!Hugo_Symbol %in% c("TTN", "MUC6", "MUC2", "MUC17", "MUC4"))
maf_f3.2_t1t0_caseinv@data = maf_f3.2_t1t0_caseinv@data %>% dplyr::filter(!Hugo_Symbol %in% c("TTN", "MUC6", "MUC2", "MUC17", "MUC4"))
maf_f3.2_t1t0_control@data = maf_f3.2_t1t0_control@data %>% dplyr::filter(!Hugo_Symbol %in% c("TTN", "MUC6", "MUC2", "MUC17", "MUC4"))

write_rds(maf_f3.2_t1t0, "maf_f3.2_t1t0.rds")
write_rds(maf_f3.2_t1t0_case, "maf_f3.2_t1t0_case.rds")
write_rds(maf_f3.2_t1t0_caseinv, "maf_f3.2_t1t0_caseinv.rds")
write_rds(maf_f3.2_t1t0_control, "maf_f3.2_t1t0_control.rds")

write_tsv(df_maf_f3.2_t1t0, "df_maf_f3.2_t1t0.maf")
write_tsv(df_maf_f3.2_t1t0_case, "df_maf_f3.2_t1t0_case.maf")
write_tsv(df_maf_f3.2_t1t0_caseinv, "df_maf_f3.2_t1t0_caseinv.maf")
write_tsv(df_maf_f3.2_t1t0_control, "df_maf_f3.2_t1t0_control.maf")


# mutsigCV -------
# module load matlab_/v81;cd ~/apps/mutsigcv/mutsig2cv;./MutSig2CV ~/projects2/ss_tnbc/testing/df_maf_f3.2_t1t0.maf ~/projects2/ss_tnbc/testing/

# Process mutsig results 
# laml.mutsig = system.file("extdata", "/mutsigcv/sig_genes.txt", package = "maftools")

ms = read.csv(file ="/home/rstudio10/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/02_annotate_mutations/mutsigcv_all/sig_genes.txt", 
              stringsAsFactors = FALSE, sep = "\t")
ms = ms[,c("gene", "q")]
#ms$q = -log10(ms$q)
ms1 = ms %>% dplyr::filter(q<0.05)
gene_list = ms1$gene
oncoplot(maf = maf_f3.2_t1t0, sortByAnnotation = "case_contrl")

genes_list = c("TP53", "PIK3CA", "DDX6", "NOTCH2","PLEC", "ERBB2", "KMT2C", "AKT1","TMEM132C", "ATM", "EGFR", "FLG2", "FOXA1",
               "TENM2", "ABCA13", "CDH23", "BTBD10", "XPO1", "DYNC1H1", "IL8", "DMGDH", "FBXW11", "SLC24A1", "MEF2A", "VHL", "SMAD4", "DDX")
gene_list2 = c(gene_list, genes_list) %>% unique()

# TMB -------
col = c("Missense_Mutation" = "purple", "Nonsense_Mutation" = "gold", "Silent" = "tomato", "Nonstop_Mutation" = "darkgreen", "Multi_Hit" = "black")
col_clin = list(c("case" = "pink", "control" = "purple"))
col_clin2 = c("case_dcis" = "coral","case_invasive" = "coral4", "control" = "seagreen2")
col_rec = c("dcis" = "mediumorchid1", "invasive" = "mediumorchid4")
fabcolors = c("orange", "purple", "grey")
names(fabcolors) = c("case", "control", "unknown")


sourcecolors = RColorBrewer::brewer.pal(n = 4,name = 'Set3')
names(sourcecolors) = c("duke", "mda", "nki", "sloane")

clin_colors = list(case_cntrl = fabcolors, source = sourcecolors, case_cntrl_new = col_clin2, recurrence = col_rec)


maf_f3.2_t1t0.tmb = maftools::tmb(maf = maf_f3.2_t1t0, logScale = FALSE)

df_tmb = maf_f3.2_t1t0.tmb
df_tmb = left_join(df_tmb, sample_idx)

#plot titv summary
pdf(file = "tmb_summary_v0.pdf", width = 9, height = 4)
maf_f3.2_t1t0_c = maftools::tmb(maf = maf_f3.2_t1t0_case, logScale = TRUE)
maf_f3.2_t1t0_cinv = maftools::tmb(maf = maf_f3.2_t1t0_caseinv, logScale = TRUE)
maf_f3.2_t1t0_co = maftools::tmb(maf = maf_f3.2_t1t0_control, logScale = TRUE)
dev.off()

pdf(file = "lollipop_v0.pdf", width = 9, height = 4)
lollipopPlot(maf_f3.2_t1t0_control, gene = "HSPG2")
lollipopPlot(maf_f3.2_t1t0_caseinv, gene = "HSPG2")
lollipopPlot(maf_f3.2_t1t0_control, gene = "TRIO")
lollipopPlot(maf_f3.2_t1t0_caseinv, gene = "APC")
lollipopPlot(maf_f3.2_t1t0_control, gene = "APC")

df_tmb$newvar = paste0(df_tmb$case_cntrl, "_", df_tmb$recurrence)
df_tmb = df_tmb %>% mutate(newvar = if_else(newvar == "control_none", "control", newvar)) 
df_tmb$newvar  = factor(df_tmb$newvar, levels = c("control", "case_dcis", "case_invasive") )
ggboxplot(data = df_tmb, x = "newvar", y = "total_perMB")
my_comparisons <- list( c("control", "case_dcis"), c("control", "case_invasive"), c("case_dcis", "case_invasive"))
p0 = ggboxplot(df_tmb, x = "newvar", y = "total_perMB",
               fill = "newvar", palette = col_clin2,
               add = "jitter") + scale_y_log10() + 
        xlab("case_cntrl") + ylab("Log10(total_perMB)") + ggtitle(paste0("Mutation Burden in Primary timepoint")) +
        theme(title = element_text(size = 8, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt"))) + 
        stat_compare_means(comparisons = my_comparisons);p0
ggsave("TMB_boxplot_case_control.png", plot = p0, width = 4, height = 6)


pdf(file = "plotmafSummary_v0.pdf", width = 15, height = 15)
plotmafSummary(maf = maf_f3.2_t1t0, rmOutlier = FALSE, addStat = 'median', dashboard = TRUE, 
               titvRaw = FALSE, showBarcodes = T, top = 20, color = col, textSize = 0.3)
dev.off()


pdf(file = "oncoplot_v0.pdf", width = 12, height = 8)
oncoplot(maf = maf_f3.2_t1t0_control, genes = gene_list2, annotationColor = clin_colors, removeNonMutated = T,
         clinicalFeatures = c("source", "case_cntrl", "recurrence"), sortByAnnotation = FALSE)
oncoplot(maf = maf_f3.2_t1t0_case, genes = gene_list2, annotationColor = clin_colors, removeNonMutated = T,
         clinicalFeatures = c("source", "case_cntrl", "recurrence"), sortByAnnotation = FALSE)
oncoplot(maf = maf_f3.2_t1t0_caseinv, genes = gene_list2, annotationColor = clin_colors, removeNonMutated = T,
         clinicalFeatures = c("source", "case_cntrl", "recurrence"), sortByAnnotation = FALSE)
dev.off()


# cohort comparison
casedcis_vs_control = mafCompare(maf_f3.2_t1t0_case, maf_f3.2_t1t0_control, m1Name = "case", m2Name = "controls")
casedcis_vs_control$results = casedcis_vs_control$results %>% dplyr::filter((!Hugo_Symbol %in% c("TTN", "MUC6", "MUC2", "MUC17", "MUC4")))

case_vs_case = mafCompare(maf_f3.2_t1t0_case, maf_f3.2_t1t0_caseinv, m1Name = "case-dcis", m2Name = "case-invasive")
case_vs_case$results = case_vs_case$results %>% dplyr::filter((!Hugo_Symbol %in% c("TTN", "MUC6", "MUC2", "MUC17", "MUC4")))

control_vs_caseinv = mafCompare(maf_f3.2_t1t0_control, maf_f3.2_t1t0_caseinv, m1Name = "control", m2Name = "case-invasive")
control_vs_caseinv$results = control_vs_caseinv$results %>% dplyr::filter((!Hugo_Symbol %in% c("TTN", "MUC6", "MUC2", "MUC17", "MUC4")))

pdf(file = "forestplot_v0.pdf", width = 5, height = 6)
forestPlot(mafCompareRes = casedcis_vs_control, fdr = 0.2)
forestPlot(mafCompareRes = case_vs_case, fdr = 0.2)
forestPlot(mafCompareRes = control_vs_caseinv, fdr = 0.2)
dev.off()

oncoplot(maf = maf_f3.2_t1t0_case, sortByAnnotation = "case_contrl")

gene_list3 = c("TP53", "PIK3CA", "APC", "SYNE1", "ATM",  "ZDHHC11",  "CDH23","PLEC","FLG2",  "AKT1", "TMEM132C",  "ERBB2", "CPAMD8",
               "KMT2C", "TENM2", "EGFR", "FOXA1", "ABCA13")
gene_list4 = c(gene_list3, gene_list2) %>% unique()

pdf(file = "co_oncoplot_v0.pdf", width = 12, height = 6)
coOncoplot(m1 = maf_f3.2_t1t0_case, m2 = maf_f3.2_t1t0_control, genes = gene_list3, m1Name = 'cases', m2Name = 'controls',removeNonMutated = T, 
           clinicalFeatures1 = "source", clinicalFeatures2 = "source", keepGeneOrder = F)
dev.off()

#somaticInteractions
pdf(file = "somatic_interactions_v0.pdf", width = 6, height = 6)
somaticInteractions(maf = maf_f3.2_t1t0_control, top = 20, pvalue = c(0.05, 0.1))
somaticInteractions(maf = maf_f3.2_t1t0_case, top = 20, pvalue = c(0.05, 0.1))
somaticInteractions(maf = maf_f3.2_t1t0_caseinv, top = 20, pvalue = c(0.05, 0.1))
dev.off()

#top vafs
pdf(file = "topvafs_v0.pdf", width = 4, height = 4)
plotVaf(maf = maf_f3.2_t1t0_control, vafCol = 'i_TumorVAF_WU', top = 12)
plotVaf(maf = maf_f3.2_t1t0_case, vafCol = 'i_TumorVAF_WU', top = 12)
plotVaf(maf = maf_f3.2_t1t0_caseinv, vafCol = 'i_TumorVAF_WU', top = 12)
dev.off()

pdf(file = "oncegenicpathways_v0.pdf", width = 4, height = 3)
OncogenicPathways(maf = maf_f3.2_t1t0_control)
OncogenicPathways(maf = maf_f3.2_t1t0_case)
OncogenicPathways(maf = maf_f3.2_t1t0_caseinv)
dev.off()

PlotOncogenicPathways(maf = maf_f3.2_t1t0_case, pathways = "RTK-RAS")
PlotOncogenicPathways(maf = maf_f3.2_t1t0_control, pathways = "RTK-RAS")

# mutational signatures
# Mutation sigs -------
install.packages("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome_casecontrol/analysis/02_annotate_mutations/BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz", repos = NULL, type="source")
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
maf_f3.2_t1t0_rem8.tnm = trinucleotideMatrix(maf = maf_f3.2_t1t0, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
maf_f3.2_t1t0_rem8.tnm_case = trinucleotideMatrix(maf = maf_f3.2_t1t0_case, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
maf_f3.2_t1t0_rem8.tnm_control = trinucleotideMatrix(maf = maf_f3.2_t1t0_control, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
# Above function performs two steps:
# Estimates APOBEC enrichment scores
# Prepares a mutational matrix for signature analysis.
library('NMF')
library(pheatmap)

maf_f3.2_t1t0_rem8.tnm_sig = estimateSignatures(mat = maf_f3.2_t1t0_rem8.tnm, nTry = 6)
maf_f3.2_t1t0_rem8.tnm_case.sig = estimateSignatures(mat = maf_f3.2_t1t0_rem8.tnm_case, nTry = 6)
maf_f3.2_t1t0_rem8.tnm_control.sig = estimateSignatures(mat = maf_f3.2_t1t0_rem8.tnm_control, nTry = 6)

plotCophenetic(res = maf_f3.2_t1t0_rem8.tnm_sig) #n=4
plotCophenetic(res = maf_f3.2_t1t0_rem8.tnm_case.sig) #n=4
plotCophenetic(res = maf_f3.2_t1t0_rem8.tnm_case.sig)#n=4

maf_f3.2_t1t0_rem8.tnm.sig = extractSignatures(mat = maf_f3.2_t1t0_rem8.tnm, n = 4)
maf_f3.2_t1t0_rem8.tnm_case.sig = extractSignatures(mat = maf_f3.2_t1t0_rem8.tnm_case, n = 4)
maf_f3.2_t1t0_rem8.tnm_control.sig = extractSignatures(mat = maf_f3.2_t1t0_rem8.tnm_control, n = 4)

pdf(file = file.path("mut_sigs_rem8_v0.pdf"), width = 9, height = 5)
plotSignatures(maf_f3.2_t1t0_rem8.tnm.sig, title_size = 0.8, contributions = TRUE, show_barcodes = TRUE, font_size = 0.6, axis_lwd = 1)
dev.off()
pdf(file = file.path("mut_sigs_rem8_v0_CASE.pdf"), width = 9, height = 5)
plotSignatures(maf_f3.2_t1t0_rem8.tnm_case.sig, title_size = 0.8, contributions = TRUE, show_barcodes = TRUE, font_size = 0.6, axis_lwd = 1)
dev.off()

pdf(file = file.path("mut_sigs_rem8_v0_CONTROL.pdf"), width = 9, height = 5)
plotSignatures(maf_f3.2_t1t0_rem8.tnm_control.sig, title_size = 0.8, contributions = TRUE, show_barcodes = TRUE, font_size = 0.6, axis_lwd = 1)
dev.off()

p = pheatmap(mat = maf_f3.2_t1t0_rem8.tnm.sig$coSineSimMat, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
library('pheatmap')
pdf(file = file.path("heatmap_mut_sigs_rem8_v1.pdf"), width = 8, height = 4)
p
dev.off()


# mut sigs deconstruc sigs ------
colors_neon = c("aquamarine1", "bisque", "blue", "brown", "brown1", "cadetblue", "cyan", "chartreuse3", "chocolate", "coral1", "darkorange", 
                "cornflowerblue", "darkgoldenrod", "darkolivegreen", "darkmagenta", "darkolivegreen1", "deeppink1", "grey", "purple", "darkblue", "black", "pink1", "plum1", "yellow1", "olivedrab1",
                "coral4", "hotpink", "tan1")
colourCount = 50
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colors_set1 = getPalette(colourCount)
colors_dark = c("blue", "red", "green", "gold2", "navy", "magenta", "darkgreen","deepskyblue",
                "orange", "deeppink", "springgreen4","darkturquoise", "brown3", "purple",
                "yellowgreen", "plum3", "darkgrey", "lavender", "cyan", "forestgreen", "indianred1", "lightblue", "slategray2", "cornsilk3", "black", colors_neon)

df_maf_v3 = df_mrg6
df_maf_v3$Tumor_Sample_Barcode = df_maf_v3$sample_idx
df_maf_v3$chr = paste0("chr", df_maf_v3$chr)

obj_mut = deconstructSigs::mut.to.sigs.input(mut.ref = as.data.frame(df_maf_v3), 
                                             sample.id = "Tumor_Sample_Barcode", 
                                             chr = "chr", 
                                             pos = "start", 
                                             ref = "ref_allele", 
                                             alt = "alt_allele", 
                                             bsg = BSgenome.Hsapiens.UCSC.hg19)

trinucs_selected <- obj_mut[rowSums(obj_mut)>50,]

#Run weight estimation

#DeconstructSigs runs on a per-sample basis, hence we have to call whichSignatures on each sample separately with the input as indicated below.

# initialize a list of the length of samples 
results <- vector("list", nrow(trinucs_selected))
names(results) <- row.names(trinucs_selected)
# run the estimation of exposures for each sample and save the results in the list
for(sID in row.names(trinucs_selected) ){
        results[[sID]] <- deconstructSigs::whichSignatures(trinucs_selected, # the matrix generated with mut.to.sigs.input 
                                                           sample.id=sID, # the current sample ID
                                                           signatures.ref=signatures.cosmic, # the data.frame with the signatures that comes with deconstructSigs
                                                           tri.counts.method="exome2genome", # which normalization method to use
                                                           contexts.needed=TRUE) # set to TRUE if your input matrix contains counts instead of frequencies
}

# convert the exposures for each sample into a sample x signatures matrix
expo <- do.call("rbind", sapply(results, "[", 1))
# add the unknown value to the matrix such that the contributions add up to 1 per sample
Signature.unknown <- unlist(sapply(results, "[", 5))
expo <- cbind(expo, Signature.unknown)
rownames(expo) = vapply(strsplit(rownames(expo), ".weights", fixed = TRUE), "[", "", 1)

# reorder samples by similarity in their signature profiles
onco_order_2 = onco_order_1[onco_order_1 %in%rownames(expo)]
# trick base graphics into putting the legend outside of the plot
par(mar=c(2,2), xpd=TRUE)
barplot(t(as.matrix(expo)), las=2, col=colors_dark)
legend("topright", inset=c(-0.1,0), legend=1:31, fill=colors_dark, title="Signature", ncol=2)

# play with df
df_res = expo
df_res$samples = rownames(df_res)
#df_res$samples = vapply(strsplit(df_res$samples, ".weights", fixed = TRUE), "[", "", 1)

# column annotation
col_anno = data.frame(pt = df_res$samples)
col_anno = left_join(col_anno, sample_idx1, by = c("pt" = "Tumor_Sample_Barcode"))

df_res = df_res[coldata1$pt,]
coldata = data.frame(col_anno, 
                     stringsAsFactors = F)
coldata1 = coldata %>% arrange(case_cntrl)
table(coldata1$pt == colnames(t(df_res))) # check sequence
ha_top = coldata1 %>% select(case_cntrl, source) %>% HeatmapAnnotation(df = ., show_annotation_name = TRUE, 
                                                                       col = clin_colors)
df_res = df_res[, -32]
sig_zero = colnames(df_res)[colSums(df_res) > 0]
df_res1 = df_res[,colnames(df_res) %in% sig_zero];dim(df_res1) #[1] 37 22
ht_v2 = Heatmap(matrix = t(df_res1), name = "weights",top_annotation = ha_top,
                row_names_gp = gpar(cex = 0.7),cluster_columns = TRUE, 
                column_names_gp = gpar(cex = 0.7),column_dend_height = unit(10, "mm"))

pdf(file  = "heatmap_mut_sigs_ordered_sigs.pdf", width = 8, height = 5)
ht_v2
dev.off() 

top_sigs = table(colnames(expo)[apply(expo, 1, which.max)])

library(ggplot2)
library(reshape)
toplot <- melt(as.matrix(expo))
names(toplot) <- c("sample", "signature", "weight")

top_sigs = table(colnames(expo)[apply(expo, 1, which.max)])


# order samples by similarity
toplot$sample <- factor(toplot$sample, levels=row.names(expo)[hclust(dist(expo))$order])
toplot = left_join(toplot, coldata, c("sample" = "pt"))

# plot
ggplot(data=toplot, aes(x=sample, y=weight, fill=signature)) + geom_bar(stat="identity") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_fill_manual(values = colors_dark)

ggplot(data=toplot, aes(x=sample, y=weight, fill=signature)) + geom_bar(stat="identity") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_fill_manual(values = colors_dark)
toplot1 = toplot %>% dplyr::filter(signature %in% c("Signature.1", "Signature.2", "Signature.3",
                                                    "Signature.5", "Signature.6", "Signature.12", "Signature.30", "Signature.13", "Signature.16")) 
g1 = ggplot(data = toplot1, aes(x = case_cntrl, y=weight, fill=signature)) + 
        #facet_wrap(~."signature") + 
        geom_bar(stat = "identity") + 
        theme_cowplot() + 
        theme(axis.text.x = element_text(angle = 90))+
        scale_fill_manual(values = colors_dark)
save_plot("barplot_case_nctrl_weights.pdf", plot = g1, base_height = 4, base_width = 3)

g2 = ggplot(data = toplot1, aes(x = source, y=weight, fill=signature)) + 
        #facet_wrap(~."signature") + 
        geom_bar(stat = "identity") + 
        theme_cowplot() + 
        theme(axis.text.x = element_text(angle = 90))+
        scale_fill_manual(values = colors_dark)
save_plot("barplot_source_weights.pdf", plot = g2, base_height = 4, base_width = 3)


toplot2 = toplot1
toplot2$signature1 = as.character(toplot2$signature)
toplot2$signature1 = vapply(strsplit(toplot2$signature1, "Signature.", fixed = TRUE), "[", "", 2) 
toplot2$signature1 = as.numeric(toplot2$signature1)

toplot2 = toplot2$sample %>% arrange(case_cntrl) 

g3 = ggplot(data=toplot2, aes(x=sample, y=weight, fill=signature)) + 
        geom_bar(stat="identity") + 
        theme_cowplot() + coord_flip() + 
        scale_fill_manual(values = colors_dark) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

g_row = ggplot(data=toplot2, aes(x=sample, fill=case_cntrl)) + 
        geom_bar(position="identity") + coord_flip() + 
        theme_cowplot() + ylab("") + xlab("") + 
        scale_fill_manual(values = c("lightblue", "pink", "grey")) + 
        theme(axis.text.y = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())

save_plot("barplot_sample_weights.pdf", plot = g3, base_height = 8, base_width = 5)
save_plot("barplot_case_cntrl_frm_mutfile_weights.pdf", plot = g_row, base_height = 8, base_width = 3)

g3 + g_row + g_row1

View(table(df_mrg6$sample_idx))

# mut count ------
my_comparisons <- list( c("control", "case_dcis"), c("control", "case_invasive"), c("case_dcis", "case_invasive"))
p0 = ggboxplot(df_tmb, x = "newvar", y = "total_perMB",
               fill = "newvar", palette = col_clin2,
               add = "jitter") + scale_y_log10() + 
        xlab("case_cntrl") + ylab("Log10(total_perMB)") + ggtitle(paste0("Mutation Burden in Primary timepoint")) +
        theme(title = element_text(size = 8, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt"))) + 
        stat_compare_means(comparisons = my_comparisons);p0
ggsave("TMB_boxplot_case_control.png", plot = p0, width = 4, height = 6)

df_mrg5$newvar = paste0(df_mrg5$case_cntrl, "_", df_mrg5$recurrence)
df_mrg5 = df_mrg5 %>% mutate(newvar = if_else(newvar == "control_none", "control", newvar)) 
df_mrg5$newvar  = factor(df_mrg5$newvar, levels = c("control", "case_dcis", "case_invasive") )

df4_box = df_mrg5 %>% dplyr::select(newvar, sample_idx) %>%  
        group_by(sample_idx) %>% mutate(num_mut_samp = n(), tmb = num_mut_samp/50) %>% ungroup %>% distinct(newvar, sample_idx, num_mut_samp)

p0 = ggboxplot(df4_box, x = "newvar", y = "num_mut_samp",
               fill = "newvar", palette = col_clin2,
               add = "jitter") + scale_y_log10() + 
        xlab("Recurrence Type") + ylab("Log10(Mutation count)") + ggtitle("Mutations in Primary timepoint") +
        theme(title = element_text(size = 8, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt"))) + 
        stat_compare_means(comparisons = my_comparisons);p0
ggsave("mutation_count.png", plot = p0, width = 4, height = 6)
df4_box_arr = df4_box %>% dplyr::arrange(newvar, num_mut_samp) 
        
p0 = ggbarplot(df4_box_arr, x = "sample_idx", y = "num_mut_samp",
               fill = "newvar", palette = col_clin2) +  
        xlab("Recurrence Type") + ylab("Log10(Mutation count)") + ggtitle("Mutations in Primary timepoint") +
        theme(title = element_text(size = 8, face = "bold"),
              axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")));p0
ggsave("mutation_count_patients.png", plot = p0, width = 10, height = 6)


p0 = ggboxplot(df_mrg5, x = "newvar", y = "tumor_f", add.params = list(size = 0.1, jitter = 0.2), 
               #order = c("private_primary", "shared_primary","private_recurrence", "shared_recurrence"),
               color = "newvar", palette = col_clin2,#fill = "patient"
               add = "jitter") + scale_y_continuous(trans='log2') +
        xlab("Mutation Type") + ylab("log2(VAF)") + ggtitle("VAF for cases vs controls mutations") +
        theme(title = element_text(size = 8, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt"))) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test");p0

ggsave("VAF.png", plot = p0, width = 4, height = 6)
df_mrg5_box = df_mrg5 %>% dplyr::select(newvar, sample_idx, tumor_f) %>%  
        group_by(sample_idx) %>% mutate(avf_tvaf = mean(tumor_f)) %>% ungroup %>% distinct(newvar, sample_idx, avf_tvaf)

p0 = ggboxplot(df_mrg5_box, x = "newvar", y = "avf_tvaf", add.params = list(size = 0.7, jitter = 0.2), 
               #order = c("private_primary", "shared_primary","private_recurrence", "shared_recurrence"),
               color = "newvar", palette = col_clin2,#fill = "patient"
               add = "jitter") + scale_y_continuous(trans='log2') +
        xlab("Mutation Type") + ylab("log2(VAF)") + ggtitle("Avg. VAF for cases vs controls patients") +
        theme(title = element_text(size = 8, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt"))) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test");p0

ggsave("VAF_avg.png", plot = p0, width = 4, height = 6)
