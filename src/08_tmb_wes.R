# Maria Roman Escorza - 2022 12 

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')
system('mkdir ./results/tmb')

library(maftools)
library(ggplot2)
library(ggpubr)


# Load WES filtered mutations ---------------------------------------------

eventDataFrame_mutect <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds')
eventDataFrame_indel <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds')

eventDataFrame_mutect <- eventDataFrame_mutect[,c('gene.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                  'exonicfunc.knowngene', 'exonicfunc.knowngene', 'patient_id', 'first_subseq_event')]

eventDataFrame_indel <- eventDataFrame_indel[,c('gene.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                  'exonicfunc.knowngene', 'exonicfunc.knowngene', 'patient_id', 'first_subseq_event')]

eventDataFrame_all <- rbind(eventDataFrame_mutect, eventDataFrame_indel)


# Prepare maf file --------------------------------------------------------

# change header
colnames(eventDataFrame_all) <- c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'first_subseq_event')

# rename  mutations
eventDataFrame_all <- eventDataFrame_all[!eventDataFrame_all$Variant_Classification %in% c("synonymous SNV"),]
eventDataFrame_all$Variant_Classification[eventDataFrame_all$Variant_Classification %in% c("frameshift deletion")] <- "Frame_Shift_Del"
eventDataFrame_all$Variant_Classification[eventDataFrame_all$Variant_Classification %in% c("frameshift insertion")] <- "Frame_Shift_Ins"
eventDataFrame_all$Variant_Classification[eventDataFrame_all$Variant_Classification %in% c("nonframeshift deletion")] <- "In_Frame_Del"
eventDataFrame_all$Variant_Classification[eventDataFrame_all$Variant_Classification %in% c("nonframeshift insertion")] <- "In_Frame_Ins"
eventDataFrame_all$Variant_Classification[eventDataFrame_all$Variant_Classification %in% c("nonsynonymous SNV")] <- "Missense_Mutation"
eventDataFrame_all$Variant_Classification[eventDataFrame_all$Variant_Classification %in% c("stoploss", "stopgain")] <- "Nonsense_Mutation"
eventDataFrame_all$Variant_Classification[eventDataFrame_all$Variant_Classification %in% c(".")] <- "Spice_Site"

eventDataFrame_all$Variant_Type[eventDataFrame_all$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Spice_Site")] <- "SNP"
eventDataFrame_all$Variant_Type[eventDataFrame_all$Variant_Classification %in% c("Frame_Shift_Ins","In_Frame_Ins")] <- "INS"
eventDataFrame_all$Variant_Type[eventDataFrame_all$Variant_Classification %in% c("Frame_Shift_Del","In_Frame_Del")] <- "DEL"

# define cases and controls
eventDataFrame_all$first_subseq_event <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_case <- eventDataFrame_all[eventDataFrame_all$first_subseq_event == 'case',-ncol(eventDataFrame_all)]
eventDataFrame_control <- eventDataFrame_all[eventDataFrame_all$first_subseq_event == 'control',-ncol(eventDataFrame_all)]

# read data frame as maf
maf_case <- read.maf(eventDataFrame_case)
maf_control <- read.maf(eventDataFrame_control)


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

