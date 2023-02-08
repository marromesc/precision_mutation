# Maria Roman Escorza - 2022 11 30 

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

library(dplyr)
library(data.table)
library(biomaRt)

source('./lib/readMetadata.R')
source('./lib/MyNumeric.R')
source('./lib/GetCosmicNumber.R')
source('./lib/filterByBed.R')

openclinica_datapath <- './data/updated_Sloane_20221114_clinicaldata_caco_genomic_samples.xlsx'
panel_datapath <- './data/Panel/DCIS_Precision_Panel_KCL'

bed <- read.delim('./data/Panel/Sloane_Covered.bed', skip = 3, sep = "\t", header = FALSE, as.is = TRUE)

# read openclinica master sheet
openclinica <- readMetadata(openclinica_datapath)
openclinica$batch <- substr(openclinica$patient_id, start = 1, stop = 7)
openclinica <- openclinica[openclinica$batch == 'PRE_SLO',]
openclinica <- openclinica[openclinica$panel == 1,]; dim(openclinica)


# Load Sloane data --------------------------------------------------------

# read kcl panel seq cases and controls 
files <- list.files(panel_datapath, pattern = '*_bcf_processed.vcf.hg19_multianno.txt', full.names = T)
samples <- gsub("\\_.*","",gsub(".*/","",files))
ifelse ( length(samples[!(samples %in% unique(openclinica$panel_id_sloane))]) == 0, 'Everything okay', stop('Missing samples in openclinica')) # check if any sequenced sample is missing in openclinica

openclinica <- openclinica[openclinica$panel!='no_normal',]  # remaining samples aren't normal paired, then they were excluded of the analysis

rownames(openclinica) <- openclinica$panel_id_sloane

openclinica <- openclinica[samples,]; dim(openclinica)
write.table(openclinica, './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_Panel_Sloane_Samples.txt', sep = '\t', quote = FALSE, row.names = FALSE)

Ind <- which(openclinica$first_subseq_event %in% c('death', 'NA', 'ipsilateral DCIS', 'ipsilateral IBC') & openclinica$surgery_final == 'BCS')

df <- lapply(files[Ind], fread)

# name mutation lists and merge samples
names(df) <- openclinica[Ind, 'patient_id']

for (i in 1:length(df)){
  df[[i]] <- as.data.frame(df[[i]])
  df[[i]]$patient_id <- names(df)[i]
}

df <- do.call(rbind, df); dim(df)

# add clinical data
df <- merge(df, openclinica); dim(df)


# Filter somatic mutations ------------------------------------------------

df <- df %>% tidyr::separate(Otherinfo13, sep=':', into = c('GT_NOR', 'AD_NOR', 'AF_NOR', 'DP_NOR', 'F1R2_NOR', 'F2R2_NOR', 'MBQ_NOR', 'MFRL_NOR', 'MMQ_NOR', 'MPOS_NOR', 'ORIGINAL_CONTIG_MISMATCH_NOR', 'SA_MAP_AF_NOR', 'SA_POST_PROB_NOR'))
df <- df %>% tidyr::separate(Otherinfo14, sep=':', into = c('GT_PDCIS', 'AD_PDCIS', 'AF_PDCIS', 'DP_PDCIS', 'F1R2_PDCIS', 'F2R2_PDCIS', 'MFRL_PDCIS', 'MMQ_PDCIS', 'MPOS_PDCIS', 'ORIGINAL_CONTIG_MISMATCH_PDCIS', 'SA_MAP_AF_PDCIS', 'SA_POST_PROB_PDCIS'))

ifelse(length(Ind) %in% length(unique(df$patient_id)), 'Everything is okay', stop("There are samples which weren't read"))
write.table(df, './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect.txt', sep = '\t', quote = FALSE, row.names = FALSE)

# filter variants in gene panel
df <- filterByBed(anned=df, bed=bed, chrom = 'Chr', start = 'Start', end = 'End')

# filter variants to have a minimum depth of 30 reads
df <- df[as.numeric(df$DP_PDCIS) >= 30,]; dim(df)

# removing synonymous snps
df <- df[!df$ExonicFunc.refGene %in% c("synonymous SNV", "unknown"),]; dim(df)

# filtering exonic and splicing
df <- df[df$Func.refGene %in% c('splicing', 'exonic', 'exonic;splicing'),]; dim(df)
df$Func.refGene <- ifelse(df$Func.refGene == 'exonic;splicing', 'exonic', df$Func.refGene) # "splicing" in ANNOVAR is defined as variant that is within 2-bp away from an exon/intron boundary by default, but the threshold can be changed by the --splicing_threshold argument. Before Feb 2013, if "exonic,splicing" is shown, it means that this is a variant within exon but close to exon/intron boundary; this behavior is due to historical reason, when a user requested that exonic variants near splicing sites be annotated with splicing as well. However, I continue to get user emails complaining about this behavior despite my best efforts to put explanation in the ANNOVAR website with details. Therefore, starting from Feb 2013 , "splicing" only refers to the 2bp in the intron that is close to an exon, and if you want to have the same behavior as before, add -exonicsplicing argument.

# filter mutations in esp6500
df <- df[MyNumeric(df$esp6500siv2_all)<0.01,]; dim(df)

# filter mutations in exac
df <- df[MyNumeric(df$ExAC_ALL)<0.01,]; dim(df)

# filter mutations in 100G
df <- df[MyNumeric(df$ALL.sites.2015_08)<0.01,]; dim(df)

# filter mutations in cosmic
df <- df[!(df$Ref=="C" & df$Alt=="T" & GetCosmicNumber(df$cosmic70)<3 & as.numeric(df$AF_PDCIS)<0.1),]; dim(df)
df <- df[!(df$Ref=="G" & df$Alt=="A" & GetCosmicNumber(df$cosmic70)<3 & as.numeric(df$AF_PDCIS)<0.1),]; dim(df)

# keep splice variants if +5 of gene start and end
grch37 = useEnsembl(biomart="genes",GRCh=37, dataset="hsapiens_gene_ensembl")
gene_annotation <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                        'start_position', 'end_position', 'band','external_synonym'),
                         filters = 'hgnc_symbol', 
                         values = unique(bed$V4), 
                         mart = grch37)
gene_annotation_2 <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                        'start_position', 'end_position', 'band','external_synonym'),
                         filters = 'external_synonym', 
                         values = unique(bed$V4), 
                         mart = grch37)
gene_annotation <- rbind(gene_annotation, gene_annotation_2)
gene_annotation <- rbind(gene_annotation, data.frame(hgnc_symbol = 'MRTFA', chromosome_name = 'chr22', start_position=40806285, end_position=41032723, band=NA, external_synonym='MKL1'))
gene_annotation$hgnc_symbol <- ifelse(gene_annotation$hgnc_symbol=='MLLT4', 'AFDN', gene_annotation$hgnc_symbol)
gene_annotation <- gene_annotation[gene_annotation$chromosome_name %in% c(1:22, 'X'),]

Ind <- which(df$Func.refGene == 'splicing')
to_remove <- c()
for(i in Ind){
  gene <- strsplit(df$Gene.refGene[i], ';')[[1]][1]
  query <- gene_annotation[gene_annotation$hgnc_symbol==gene,]
  if(nrow(query)==0){
    query <- gene_annotation[gene_annotation$external_synonym==gene,]
  } else if (nrow(query) > 1){
    query <- query[1,]
  }
  if (!isTRUE(df$Start[i] > query$start_position+5 & df$End[i] < query$end_position-5)){
    to_remove <- c(to_remove,i)
  }
}

if(!is.null(to_remove)){
  df <- df[-to_remove,]
}


# variants found in somatic clinvar were included
#df <- df[!(df$ClinVarSomatic == 0 & df$ClinVarGermline != 0),]

# filter mutations with low af
df <- df[as.numeric(df$AF_PDCIS) >= 0.05,]; dim(df)

# export filtered mutations
saveRDS(df, './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect_Filtered.rds')
write.table(df, './data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect_Filtered.txt', sep = '\t', quote = FALSE, row.names = FALSE)
