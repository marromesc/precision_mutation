# Maria Roman Escorza - 2022 11 30 

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

library(dplyr)
library(readr)
library(liftOver)

source('./lib/readMetadata.R')
source('./lib/MyNumeric.R')
source('./lib/GetCosmicNumber.R')

openclinica_datapath <- './data/updated_20221114_clinicaldata_caco_genomic_samples.xlsx'
panel_datapath <- './data/Panel/DCIS_Precision_Panel_NKI/Compiled_mutations_nki.xlsx'

GenePanel <- readRDS('./data/GenesPanel.RDS')

# read openclinica master sheet
openclinica <- readMetadata(openclinica_datapath)
openclinica$batch <- substr(openclinica$patient_id, start = 1, stop = 7)
openclinica <- openclinica[openclinica$batch == 'PRE_NKI',]
openclinica <- openclinica[openclinica$panel == 1,]
write.table(openclinica, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# Load NKI data -----------------------------------------------------------

df <- as.data.frame(readxl::read_xlsx(panel_datapath))

openclinica <- openclinica[openclinica$panel == 1,]; dim(openclinica)
df <- df[df$sample_ID_DCIS %in% openclinica$panel_id,]; dim(df)

df <- merge(df, openclinica, by.x = 'sample_ID_DCIS', by.y = 'panel_id'); dim(df)

# remove IBC mutations
df <- df[df$vaf_DCIS_DNA != 'NA',]; dim(df)

#liftover
df.gr <- df[,c('CHROM', 'POS', 'REF', 'ALT')]
df.gr$CHROM <- paste0('chr', df$CHROM)
df.gr <- makeGRangesFromDataFrame(df.gr, 
                                  keep.extra.columns = T,
                                  seqnames.field="CHROM",
                                  start.field ='POS', end.field='POS'
                                  )

ch = import.chain('/mnt/albyn/common/chains/hg38ToHg19.over.chain')
df.gr = as.data.frame(liftOver(df.gr, ch))

df$POS <- df.gr$start

for(i in 1:nrow(df)){
  alt <- df$ALT[i]
  
  if(nchar(alt) > 1 & substr(alt, 1, 1) == '-'){ #deletions
    df$ALT[i] <- '-'
    df$REF[i] <- gsub("[^A-Z]","",df$ALT)
  } else if (nchar(alt) > 1 & substr(alt, 1, 1) == '+'){ #insertions
    df$REF[i] <- '-'
    df$ALT[i] <- gsub("[^A-Z]","",df$ALT)
  }
}


# Annotate mutations ------------------------------------------------------

write.table(df[,c('CHROM', 'POS', 'POS', 'REF', 'ALT')], './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect.avinput', col.names = F, row.names = F, quote = F)
system('perl /mnt/albyn/common/annovar/table_annovar.pl --buildver hg19 ./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect.avinput -protocol refGene,exac03,cosmic70,esp6500siv2_all,ALL.sites.2015_08 -operation g,f,f,f,f -polish /mnt/albyn/common/annovar/humandb')

df.new <- data.table::fread('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect.avinput.hg19_multianno.txt')
df <- cbind(df, df.new); dim(df)
write.table(df, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# Filter somatic mutations ------------------------------------------------

# filter variants in gene panel
df <- df[df$SYMBOL %in% GenePanel,]; dim(df)

# filter low quality variants
df <- df[(MyNumeric(df$vaf_DCIS_DNA)/100) > 0.1 & as.numeric(df$cov_DCIS_DNA) > 100,]; dim(df)

# removing synonymous snps
df <- df[!df$ExonicFunc.refGene %in% c("synonymous SNV"),]; dim(df)

# filtering exonic and splicing
df <- df[df$Func.refGene %in% c('splicing', 'exonic', 'exonic;splicing', ''),]; dim(df)
df$Func.refGene <- ifelse(df$Func.refGene == 'exonic;splicing', 'exonic', df$Func.refGene) # "splicing" in ANNOVAR is defined as variant that is within 2-bp away from an exon/intron boundary by default, but the threshold can be changed by the --splicing_threshold argument. Before Feb 2013, if "exonic,splicing" is shown, it means that this is a variant within exon but close to exon/intron boundary; this behavior is due to historical reason, when a user requested that exonic variants near splicing sites be annotated with splicing as well. However, I continue to get user emails complaining about this behavior despite my best efforts to put explanation in the ANNOVAR website with details. Therefore, starting from Feb 2013 , "splicing" only refers to the 2bp in the intron that is close to an exon, and if you want to have the same behavior as before, add -exonicsplicing argument.

# variants with more than 5 entries in gnomad were removed
df <- df[MyNumeric(df$GNOMAD_AC)<5,]; dim(df)

# variants with more than 5 entries in gonl were removed
df <- df[MyNumeric(df$GONL_AC)<5,]; dim(df)

# filter mutations in esp6500
df <- df[MyNumeric(df$esp6500)<0.01,]; dim(df)

# filter mutations in exac
df <- df[MyNumeric(df$ExAC_ALL)<0.01,]; dim(df)

# filter mutations in 100G
df <- df[MyNumeric(df$ALL.sites.2015_08)<0.01,]; dim(df)

# filter mutations in cosmic
df <- df[!(df$REF=="C" & df$ALT=="T" & GetCosmicNumber(df$cosmic70)<3 & (as.numeric(df$vaf_DCIS_DNA)/100)<0.1),]
df <- df[!(df$REF=="G" & df$ALT=="A" & GetCosmicNumber(df$cosmic70)<3 & (as.numeric(df$vaf_DCIS_DNA)/100)<0.1),]
dim(df)

# variants found in somatic clinvar were included
df <- df[!(df$ClinVarSomatic == 0 & df$ClinVarGermline != 0),]; dim(df)

# export filtered mutations
write.table(df, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)
saveRDS(df, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.rds')


