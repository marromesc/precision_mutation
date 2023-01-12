# Maria Roman Escorza - 2022 12 

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')
system('mkdir ./results/tmb')

library(maftools)
library(ggplot2)
library(ggpubr)
library(readr)

maf_datapath <- './results/Filtered_Mutations_Compiled.maf'
meta_datapath <- './results/SampleSheet.csv'

SampleSheet <- read.csv(meta_datapath)

mafDF <- as.data.frame(read_tsv(maf_datapath))
mafDF <- mafDF[mafDF$Center %in% c('NKI1', 'NKI2', 'NKI4', 'NKI3', 'SLO1', 'SLO2', 'SLO3', 'Duke1'),]

maf_case <- read.maf(mafDF[mafDF$Entrez_Gene_Id %in% SampleSheet$patient_id[SampleSheet$case_control=='case'],])
maf_control <- read.maf(mafDF[mafDF$Entrez_Gene_Id %in% SampleSheet$patient_id[SampleSheet$case_control=='control'],])


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

