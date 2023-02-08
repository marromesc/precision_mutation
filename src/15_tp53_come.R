# Maria Roman Escorza - 2023 02 08

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

system('mkdir ./results/cn_muts')

library(readr)
library(ComplexHeatmap)
library(circlize)

source('./lib/somaticInteractions.R')
source('./lib/addCN2Muts.R')
source('./lib/mutCountMatrix.R')
source('./lib/oncoPlotDetails.R')

mutation_datapath <- './results/Filtered_Mutations_Compiled_5%.csv'
meta_datapath <- './results/SampleSheet.csv'
meta_cn_datapath <- '/home/maria/albyn/copynumber/precision_copynumber/results/preprocess_cn/SamplesInfo.csv'


# Load data ---------------------------------------------------------------

eventDataFrame <- read.csv(mutation_datapath)
eventDataFrame <- eventDataFrame[!is.na(eventDataFrame$case_control),]
SampleSheet <- read.csv(meta_datapath)
rownames(SampleSheet) <- SampleSheet$patient_id
SampleSheet_CN <- as.data.frame(read_tsv(meta_cn_datapath))
SampleSheet_CN <- SampleSheet_CN[SampleSheet_CN$first_subseq_event %in% c('ipsilateral IBC', 'NoData') & SampleSheet_CN$surgery_final == 'BCS',]
rownames(SampleSheet_CN) <- SampleSheet_CN$patient_id

SampleSheet<-SampleSheet[SampleSheet$patient_id%in%SampleSheet_CN$patient_id,]
SampleSheet_CN<-SampleSheet_CN[SampleSheet_CN$patient_id%in%SampleSheet$patient_id,]
SampleSheet <- SampleSheet[rownames(SampleSheet_CN),]


# Set up copy number data + mutations --------------------------------------

geneMatrix <- mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, 'TP53', annotation = T)
df <- data.frame(loh=ifelse(SampleSheet_CN$TP53_loh==1, 'loh', SampleSheet_CN$TP53_loh), 
                 mutation=mutCountMatrix(eventDataFrame = eventDataFrame, patients = SampleSheet$patient_id, 'TP53', annotation = T), 
                 caco = SampleSheet$case_control,
                 her2 = SampleSheet$her2,
                 er = SampleSheet$er,
                 grade = SampleSheet$grade)
df <- df[!is.na(df$loh),]
df$loh <- ifelse(df$loh==0, NA, df$loh)
df$TP53 <- ifelse(df$TP53==0, NA, df$TP53)


# Oncoplot ----------------------------------------------------------------

pdf(paste0("./results/cn_muts/oncoPrint_TP53_loh.pdf"), height = 5)
print(oncoPrint(t(as.matrix(df[1:2])), alter_fun = alter_fun,
                bottom_annotation = HeatmapAnnotation(df = data.frame(CaseControl=df$caco),
                                                      col=list(CaseControl=c('case'='green', 'control'='pink')))))
dev.off()


# Compare -----------------------------------------------------------------

df$co_occur <- ifelse(!is.na(df$loh) & !is.na(df$TP53), 1, 0)
write.csv(df, "./results/cn_muts/TP53_loh.csv", row.names = F)
fisher.test(table(df$caco, df$co_occur))
fisher.test(table(df$her2, df$co_occur))
fisher.test(table(df$er, df$co_occur))
fisher.test(table(df$grade, df$co_occur))

fisher.test(table(df$caco[df$er==0], df$co_occur[df$er==0]))
