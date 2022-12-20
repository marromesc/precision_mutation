library(dplyr)
library(readr)

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
    if (z==".")
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


# Load Sloane data --------------------------------------------------------

# read nki panel seq cases and controls 
files <- list.files('./data/Panel/DCIS_Precision_Panel_NKI', pattern = '*.vcf.hg19_multianno.txt', full.names = T)
samples <- gsub(".vcf.hg19_multianno.txt*","",gsub(".*/","",files))

#openclinica <- openclinica[openclinica$panel_id %in% samples,]; dim(openclinica)
case_control <- openclinica[openclinica$first_subseq_event %in% c('ipsilateral IBC', 'death', 'NA') & openclinica$panel == 1,'patient_id']
filt_mut <- data.table::fread('/mnt/albyn/maria/targeted_seq_nki/DCIS_variants_f_cleanNsimple_noCommonBRCA_precision_ids_CV.txt')
case_control_file <- unique(filt_mut$site_accessio)
filt_mut1 <- readxl::read_xlsx('/mnt/albyn/maria/copynumber/mutation_amp/data/Compiled_mutations_nki.xlsx')
case_control_file2 <- unique(filt_mut[['sample_ID_DCIS']])
    
case_control_file2 == case_control_file

case_control[!(case_control %in% case_control_file)]
case_control_file[!(case_control_file %in% case_control)]

case_control[!(case_control %in% samples)]

case_control[!(case_control %in% case_control_file)]
case_control_file[!(case_control_file %in% case_control)]


case_control[!(case_control %in% case_control_file2)]
case_control_file2[!(case_control_file2 %in% case_control)]

Ind <- which(samples %in% rownames(openclinica))

mafMat <- lapply(files[Ind], read_tsv, col_names = TRUE)

# name mutation lists and merge samples
names(mafMat) <- openclinica$patient_id

for (i in 1:length(mafMat)){
  mafMat[[i]] <- as.data.frame(mafMat[[i]])
  mafMat[[i]]$ID <- names(mafMat)[i]
}

mafMat <- do.call(rbind, mafMat); dim(mafMat)


# Filter somatic mutations ------------------------------------------------

mafMat <- mafMat %>% tidyr::separate(Otherinfo13, sep=':', into = c('GT', 'GQ', 'DP', 'FDP', 'RO', 'FRO', 'AO', 'FAO', 'AF', 'SAR', 'SAF', 'SRF', 'SRR', 'FSAR', 'FSAF', 'FSRF', 'FSRR'))

# filter variants in gene panel
GenePanel <- c("ARID1A","GATA3","PTEN","ATM","KMT2D","RB1","AKT1","CDH1","TP53","MAP2K4","NCOR1","NF1","ERBB2","BRCA1","RUNX1","CHEK2","BAP1","PIK3CA","FBXW7","MAP3K1","PIK3R1","KMT2C","NOTCH1","SF3B1","PBRM1","PDGFRA","CCND3","ESR1","ARID1B","EGFR","BRAF","FGFR1","MYC","CDKN2A","FGFR2","CCND1","ERBB3","MDM2","TBX3","BRCA2","IGF1R","CBFB","SMAD4","STK11","CCNE1") 
mafMat <- mafMat[mafMat$Gene.refGene %in% GenePanel,]; dim(mafMat)

# filter low quality variants
mafMat1 <- mafMat[!(as.numeric(mafMat$AF) < 0.1 & as.numeric(mafMat$DP) < 100 & as.numeric(mafMat$Otherinfo9 < 1000)),]; dim(mafMat1)
mafMat2 <- mafMat[as.numeric(mafMat$AF) > 0.1 & as.numeric(mafMat$DP) > 100 & as.numeric(mafMat$Otherinfo9 > 1000),]; dim(mafMat2)

# removing synonymous snps
mafMat <- mafMat[!mafMat$ExonicFunc.refGene %in% c("synonymous SNV"),]; dim(mafMat)

# filtering exonic and splicing
mafMat <- mafMat[mafMat$Func.refGene %in% c('splicing', 'exonic', 'exonic;splicing'),]; dim(mafMat)

# variants not found in gnomad were included
grepl('GNOMAD', mafMat$Otherinfo11)

mafMat$Otherinfo11

# variants not found gonl were included
////////////////

# filter mutations in esp6500
mafMat <- mafMat[MyNumeric(mafMat$esp6500si_all)<0.01,]

# filter mutations in exac
mafMat <- mafMat[MyNumeric(mafMat$ExAC_ALL)<0.01,]

# filter mutations in 100G
mafMat <- mafMat[MyNumeric(mafMat$x1kg2015aug_max)<0.01,]

# filter mutations in cosmic
mafMat <- mafMat[!(mafMat$ref=="C" & mafMat$alt=="T" & GetCosmicNumber(mafMat$cosmic70)<3 & mafMat$AF<0.1),]
mafMat <- mafMat[!(mafMat$ref=="G" & mafMat$alt=="A" & GetCosmicNumber(mafMat$cosmic70)<3 & mafMat$AF<0.1),]

# variants found in somatic clinvar were included
#mafMat <- mafMat[!(mafMat$ClinVarSomatic == 0 & mafMat$ClinVarGermline != 0),]

# filter mutations with low af
mafMat <- mafMat[mafMat$AF >= 0.05,]

# rename mutations
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c(".")] <- "splicing"
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c("nonsynonymous_SNV")] <- "missense"
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c("nonframeshift_deletion","nonframeshift_insertion","nonframeshift_substitution")] <- "inframe_indel"
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c("frameshift_deletion","frameshift_insertion","frameshift_substitution")] <- "frameshift"
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"

# export filtered mutations
write.table(mafMat, 'DCIS_Precision_Panel_Sloane_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)



