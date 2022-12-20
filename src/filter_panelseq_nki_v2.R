library(dplyr)
library(readr)
library(liftOver)

GetCosmicNumber <- function(x)
{
  x <- sapply(x,function(z){
    if (is.na(z))
    {
      return(0)
    }
    else if (z==".")
    {
      return(0)
    }
    else
    {
      z <- strsplit(z,"=")[[1]][3]
      z <- strsplit(z,",")[[1]]
      sum <- 0
      for (i in 1:length(z))
      {
        foo <- strsplit(z[i],"\\(")[[1]]
        sum <- sum + as.numeric(foo[1])
      }
      return(sum)
    }
  }) 
  return(as.numeric(x))
}

MyNumeric <- function(x)
{
  x <- sapply(x,function(z){
    if ((z=="." | is.na(z) | z == 'NA'))
    {
      return(0)
    }
    else
    {
      return(as.numeric(z))
    }
  }) 
  return(as.numeric(x))
}

setwd('/mnt/albyn/maria/precision_mutation')

# read openclinica master sheet
openclinica <- as.data.frame(readxl::read_xlsx('./data/updated_20221114_clinicaldata_caco_genomic_samples.xlsx')); dim(openclinica)
openclinica2 <- as.data.frame(readxl::read_xlsx('./data/updated_20221114_clinicaldata_caco_genomic_samples.xlsx', sheet = 2))
openclinica2 <- openclinica2[openclinica2$redcap_event == 'PRI',]
openclinica <- as.data.frame(merge(openclinica, openclinica2)); dim(openclinica)
openclinica$batch <- substr(openclinica$patient_id, start = 1, stop = 7)
openclinica <- openclinica[openclinica$batch == 'PRE_NKI',]
openclinica <- openclinica[openclinica$panel == 1,]
write.table(openclinica, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# Load NKI data --------------------------------------------------------

df <- as.data.frame(readxl::read_xlsx('./data/Panel/DCIS_Precision_Panel_NKI/Compiled_mutations_nki.xlsx'))
dim(df)

openclinica <- openclinica[openclinica$first_subseq_event %in% c('ipsilateral DCIS', 'ipsilateral IBC', 'death', 'NA') & openclinica$panel == 1,]; dim(openclinica)
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
system('perl /mnt/albyn/common/annovar/table_annovar.pl --buildver hg19 ./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect.avinput -protocol refGene,exac03,cosmic70,avsnp147,dbnsfp30a,esp6500siv2_all,ALL.sites.2015_08,gnomad_exome -operation g,f,f,f,f,f,f,f -polish /mnt/albyn/common/annovar/humandb')

df.new <- data.table::fread('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect.avinput.hg19_multianno.txt')
df <- cbind(df, df.new); dim(df)
write.table(df, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# Filter somatic mutations ------------------------------------------------

# filter variants in gene panel
GenePanel <- readRDS('./data/GenesPanel.RDS')
df <- df[df$SYMBOL %in% GenePanel,]; dim(df)

# filter low quality variants
df2 <- df[(as.numeric(df$vaf_DCIS_DNA)/100) > 0.1 & as.numeric(df$cov_DCIS_DNA) > 100,]; dim(df2)

# removing synonymous snps
df3 <- df2[!df2$ExonicFunc.refGene %in% c("synonymous SNV"),]; dim(df3)

# filtering exonic and splicing
df3 <- df3[df3$Func.refGene %in% c('splicing', 'exonic', 'exonic;splicing', ''),]; dim(df3)

# variants not found in gnomad were included
df4 <- df3[MyNumeric(df3$gnomAD_exome_ALL)<0.01,]; dim(df4)

# variants not found gonl were included
df4 <- df3[MyNumeric(df3$GONL_AF)<0.01,]; dim(df4)

# filter mutations in esp6500
df4 <- df4[MyNumeric(df4$esp6500)<0.01,]; dim(df4)

# filter mutations in exac
df5 <- df4[MyNumeric(df4$ExAC_ALL)<0.01,]; dim(df5)

# filter mutations in 100G
df6 <- df5[MyNumeric(df5$ALL.sites.2015_08)<0.01,]; dim(df6)

# filter mutations in cosmic
df7 <- df6[!(df6$REF=="C" & df6$ALT=="T" & GetCosmicNumber(df6$cosmic70)<3 & (as.numeric(df6$vaf_DCIS_DNA)/100)<0.1),]
df7 <- df6[!(df6$REF=="G" & df6$ALT=="A" & GetCosmicNumber(df6$cosmic70)<3 & (as.numeric(df6$vaf_DCIS_DNA)/100)<0.1),]
dim(df7)

# variants found in somatic clinvar were included
df8 <- df7[!(df7$ClinVarSomatic == 0 & df7$ClinVarGermline != 0),]; dim(df8)

# filter mutations with low af
df9 <- df8[(as.numeric(df8$vaf_DCIS_DNA)/100) >= 0.05,]; dim(df9)

# rename mutations
df9$Consequence[df9$Consequence %in% c(".", 'splice_acceptor_variant', 'splice_donor_variant')] <- "splicing"
df9$Consequence[df9$Consequence %in% c("nonsynonymous SNV", "missense_variant", 'missense_variant&splice_region_variant')] <- "missense"
  df9$Consequence[df9$Consequence %in% c("nonframeshift_deletion","nonframeshift_insertion","nonframeshift_substitution", 'inframe_insertion','inframe_deletion')] <- "inframe_indel"
df9$Consequence[df9$Consequence %in% c("frameshift deletion","frameshift insertion","frameshift substitution", "frameshift_variant")] <- "frameshift"
df9$Consequence[df9$Consequence %in% c("stopgain","stoploss","startgain","startloss", "stop_gained")] <- "nonsense"

# export filtered mutations
write.table(df9, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)
saveRDS(df9, './data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.rds')


