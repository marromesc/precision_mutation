library(dplyr)
library(data.table)

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
openclinica <- openclinica[openclinica$batch == 'PRE_SLO',]
openclinica <- openclinica[openclinica$panel == 1,]; dim(openclinica)

# map sloane and precision ids
mapping_file <- readxl::read_xlsx('/mnt/albyn/maria/master/Sloane-project-samples-MASTER-sheet-May-2022.xlsx', sheet = 2)
mapping_file <- as.data.frame(mapping_file[!(mapping_file$Tissue %in% c('Normal', 'DCIS recurrence', 'Inv recurrence', 'Inv recurrence Tubular', 'Inv recurrence NST', 'Inv recurrence 2')) 
                               & !is.na(mapping_file$'ng SureSelect'),]); dim(mapping_file)
ifelse ( length(unique(mapping_file$'DNA-Short')[!(unique(mapping_file$'DNA-Short') %in% unique(openclinica$panel_id))]) == 0, 'Everything okay', stop('Missing samples in openclinica')) # check if any sequenced sample is missing in openclinica
openclinica <- merge(openclinica, mapping_file[,c('Precision ID', 'Sloane ID')], by.x = 'patient_id', by.y = 'Precision ID'); dim(openclinica)
rownames(openclinica) <- openclinica$'Sloane ID'


# Load Sloane data --------------------------------------------------------

# read kcl panel seq cases and controls 
files <- list.files('./data/Panel/DCIS_Precision_Panel_KCL', pattern = '*_bcf_processed.vcf.hg19_multianno.txt', full.names = T)
samples <- gsub("\\_.*","",gsub(".*/","",files))
ifelse ( length(samples[!(samples %in% unique(openclinica$'Sloane ID'))]) == 0, 'Everything okay', stop('Missing samples in openclinica')) # check if any sequenced sample is missing in openclinica

openclinica$`Sloane ID`[!(openclinica$`Sloane ID` %in% samples)] # remaining samples aren't normal paired, then they were excluded of the analysis

openclinica <- openclinica[samples,]; dim(openclinica)
write.table(openclinica, './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_Panel_Sloane_Samples.txt', sep = '\t', quote = FALSE, row.names = FALSE)

Ind <- which(openclinica$first_subseq_event %in% c('death', 'NA', 'ipsilateral DCIS', 'ipsilateral IBC'))

mafMat <- lapply(files[Ind], fread)

# name mutation lists and merge samples
names(mafMat) <- openclinica[Ind, 'patient_id']

for (i in 1:length(mafMat)){
  mafMat[[i]] <- as.data.frame(mafMat[[i]])
  mafMat[[i]]$patient_id <- names(mafMat)[i]
}

mafMat <- do.call(rbind, mafMat); dim(mafMat)

# add clinical data
mafMat <- merge(mafMat, openclinica); dim(mafMat)

# Filter somatic mutations ------------------------------------------------

mafMat <- mafMat %>% tidyr::separate(Otherinfo13, sep=':', into = c('GT_NOR', 'AD_NOR', 'AF_NOR', 'DP_NOR', 'F1R2_NOR', 'F2R2_NOR', 'MBQ_NOR', 'MFRL_NOR', 'MMQ_NOR', 'MPOS_NOR', 'ORIGINAL_CONTIG_MISMATCH_NOR', 'SA_MAP_AF_NOR', 'SA_POST_PROB_NOR'))
mafMat <- mafMat %>% tidyr::separate(Otherinfo14, sep=':', into = c('GT_PDCIS', 'AD_PDCIS', 'AF_PDCIS', 'DP_PDCIS', 'F1R2_PDCIS', 'F2R2_PDCIS', 'MFRL_PDCIS', 'MMQ_PDCIS', 'MPOS_PDCIS', 'ORIGINAL_CONTIG_MISMATCH_PDCIS', 'SA_MAP_AF_PDCIS', 'SA_POST_PROB_PDCIS'))

ifelse(length(Ind) %in% length(unique(mafMat$patient_id)), 'Everything is okay', stop("There are samples which weren't read"))
write.table(mafMat, './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# filter variants in gene panel
bed <- read.delim('./data/Panel/Sloane_Covered.bed', skip = 3, sep = "\t", header = FALSE, as.is = TRUE)
GenePanel <- unique(bed[,4])
mafMat <- mafMat[mafMat$Gene.refGene %in% GenePanel,]; dim(mafMat)

# filter variants to have a minimum depth of 30 reads
mafMat <- mafMat[as.numeric(mafMat$DP_PDCIS) >= 30,]; dim(mafMat)

# removing synonymous snps
mafMat <- mafMat[!mafMat$ExonicFunc.refGene %in% c("synonymous SNV", "unknown"),]; dim(mafMat)

# filtering exonic and splicing
mafMat <- mafMat[mafMat$Func.refGene %in% c('splicing', 'exonic', 'exonic;splicing'),]; dim(mafMat)

# filter mutations in esp6500
mafMat <- mafMat[MyNumeric(mafMat$esp6500siv2_all)<0.01,]; dim(mafMat)

# filter mutations in exac
mafMat <- mafMat[MyNumeric(mafMat$ExAC_ALL)<0.01,]; dim(mafMat)

# filter mutations in 100G
mafMat <- mafMat[MyNumeric(mafMat$ALL.sites.2015_08)<0.01,]; dim(mafMat)

# filter mutations in cosmic
mafMat <- mafMat[!(mafMat$Ref=="C" & mafMat$Alt=="T" & GetCosmicNumber(mafMat$cosmic70)<3 & as.numeric(mafMat$AF_PDCIS)<0.1),]; dim(mafMat)
mafMat <- mafMat[!(mafMat$Ref=="G" & mafMat$Alt=="A" & GetCosmicNumber(mafMat$cosmic70)<3 & as.numeric(mafMat$AF_PDCIS)<0.1),]; dim(mafMat)

# variants found in somatic clinvar were included
#mafMat <- mafMat[!(mafMat$ClinVarSomatic == 0 & mafMat$ClinVarGermline != 0),]

# filter mutations with low af
mafMat <- mafMat[as.numeric(mafMat$AF_PDCIS) >= 0.05,]; dim(mafMat)

# rename mutations
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c(".")] <- "splicing"
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c("nonsynonymous SNV")] <- "missense"
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution")] <- "inframe_indel"
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c("frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift" 
mafMat$ExonicFunc.refGene[mafMat$ExonicFunc.refGene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"

ifelse(length(Ind) %in% length(unique(mafMat$patient_id)), 'Everything is okay', "There are samples with no mutations after filtering, as expected")

# export filtered mutations
saveRDS(mafMat, './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect_Filtered.rds')
write.table(mafMat, './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)



