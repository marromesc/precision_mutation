# Maria Roman Escorza - 2023 01 27

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')
system('mkdir ./results/mut_signatures')

library(maftools)
library(data.table)
library(dplyr)
library("BSgenome.Hsapiens.UCSC.hg19")

wes_indel_datapath <- './data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered_discovery_5%.rds'
wes_mutect_datapath <- './data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered_discovery_5%.rds'
wes_meta_datapath <- './results/SampleSheet.csv'


# Load data ---------------------------------------------------------------

SampleSheet <- as.data.frame(fread(wes_meta_datapath))
SampleSheet <- SampleSheet[SampleSheet$platform=='WES',]

eventDataFrame_mutect <- readRDS(wes_mutect_datapath) %>% dplyr::mutate(Hugo_Symbol=gene.knowngene, Entrez_Gene_Id=patient_id, Center=paste0(batch,'-WES'), NCBI_Build='hg19', Chromosome=chr, Start_Position=start, 
                                                                        End_Position=end, Strand='+', 
                                                                        Variant_Classification=ifelse(exonicfunc.knowngene=='nonsynonymous SNV', 'Missense_Mutation',
                                                                                                      ifelse(exonicfunc.knowngene%in%c('stopgain','stoploss'), 'Nonsense_Mutation',
                                                                                                             ifelse(exonicfunc.knowngene=='.', 'Splice_Site','NoData'))), 
                                                                        Variant_Type='SNP', Reference_Allele=ref_allele, Tumor_Seq_Allele1=ref_allele, Tumor_Seq_Allele2=alt_allele, dbSNP_RS=dbsnp_site,
                                                                        Tumor_Sample_Barcode=sample_name_pdcis, Matched_Norm_Sample_Barcode=sample_name_normal, 
                                                                        t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count, t_vaf=tumor_f, n_vaf=normal_af,
                                                                        SIFT=sift_pred, PolyPhen=polyphen2_hvar_pred, GMAF=x1kg2015aug_max, CLIN_SIG=clinsig, ExAC_AF=exac_all,
                                                                        COSMIC=cosmic70, ESP6500=esp6500siv2_all) %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_Position, 
                                                                                                                                    End_Position, Strand, Variant_Classification, Variant_Type, 
                                                                                                                                    Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, Tumor_Sample_Barcode, 
                                                                                                                                    Matched_Norm_Sample_Barcode, t_depth, t_ref_count, t_alt_count, n_depth, 
                                                                                                                                    n_ref_count, n_alt_count, t_vaf, n_vaf, SIFT, PolyPhen, GMAF, CLIN_SIG, 
                                                                                                                                    ExAC_AF, COSMIC, ESP6500, age_diagnosis, radiotherapy, time_at_risk,
                                                                                                                                    first_subseq_event, case_control, er, pr, her2, grade, rna_id, cnv_id, wes_id)

eventDataFrame_indel <- readRDS(wes_indel_datapath) %>% dplyr::mutate(Hugo_Symbol=gene.knowngene, Entrez_Gene_Id=patient_id, Center=paste0(batch,'-WES'), NCBI_Build='hg19', Chromosome=chr, Start_Position=start, 
                                                                      End_Position=end, Strand='+', 
                                                                      Variant_Classification=ifelse(exonicfunc.knowngene=='frameshift deletion', 'Frame_Shift_Del',
                                                                                                    ifelse(exonicfunc.knowngene%in%c('stopgain','stoploss'), 'Nonsense_Mutation',
                                                                                                           ifelse(exonicfunc.knowngene=='.', 'Splice_Site',
                                                                                                                  ifelse(exonicfunc.knowngene=='nonframeshift deletion', 'In_Frame_Del',
                                                                                                                         ifelse(exonicfunc.knowngene=='nonframeshift insertion', 'In_Frame_Ins',
                                                                                                                                ifelse(exonicfunc.knowngene=='frameshift insertion', 'Frame_Shift_Ins', 'NoData')))))), 
                                                                      Variant_Type=ifelse(ref_allele=='0' & alt_allele=='-', 'DEL',
                                                                                          ifelse(ref_allele=='-' & alt_allele!='-', 'INS', 'NoData')), Reference_Allele=ref_allele, Tumor_Seq_Allele1=ref_allele, Tumor_Seq_Allele2=alt_allele, dbSNP_RS='NoData',
                                                                      Tumor_Sample_Barcode=sample_name_pdcis, Matched_Norm_Sample_Barcode=sample_name_normal, 
                                                                      t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count, t_vaf=t_vaf, n_vaf=normal_af,
                                                                      SIFT=sift_pred, PolyPhen=polyphen2_hvar_pred, GMAF=x1kg2015aug_max, CLIN_SIG=clinsig, ExAC_AF=exac_all,
                                                                      COSMIC=cosmic70, ESP6500=esp6500siv2_all) %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_Position, 
                                                                                                                                  End_Position, Strand, Variant_Classification, Variant_Type, 
                                                                                                                                  Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, Tumor_Sample_Barcode, 
                                                                                                                                  Matched_Norm_Sample_Barcode, t_depth, t_ref_count, t_alt_count, n_depth, 
                                                                                                                                  n_ref_count, n_alt_count, t_vaf, n_vaf, SIFT, PolyPhen, GMAF, CLIN_SIG, 
                                                                                                                                  ExAC_AF, COSMIC, ESP6500, age_diagnosis, radiotherapy, time_at_risk,
                                                                                                                                  first_subseq_event, case_control, er, pr, her2, grade, rna_id, cnv_id, wes_id)

eventDataFrame_wes <- rbind(eventDataFrame_mutect, eventDataFrame_indel)

maf <- read.maf(eventDataFrame_wes)


# Estimate signatures -----------------------------------------------------

tnm = trinucleotideMatrix(maf = maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
write.csv(tnm[[2]], './results/mut_signatures/mut_signatures.csv', quote = F, row.names = F)


# Load data for Signal ----------------------------------------------------

# preparing input for signal mutational analysis https://signal.mutationalsignatures.com

eventDataFrame_mutect <- readRDS(wes_mutect_datapath) %>% dplyr::mutate('Sample name'=patient_id, Chromosome=chr, Position=start, 
                                                                        'Original base'=ref_allele, 'Mutated base'=alt_allele) %>% dplyr::select('Sample name', Chromosome, Position, 'Original base', 'Mutated base')

eventDataFrame_indel <- readRDS(wes_indel_datapath) %>% dplyr::mutate('Sample name'=patient_id, Chromosome=chr, Position=start, 
                                                                      'Original base'=ref_allele, 'Mutated base'=alt_allele) %>% dplyr::select('Sample name', Chromosome, Position, 'Original base', 'Mutated base')

eventDataFrame_wes <- rbind(eventDataFrame_mutect, eventDataFrame_indel)
eventDataFrame_wes <- eventDataFrame_wes[eventDataFrame_wes$`Sample name` == 'PRE_SLO4673',]

write.csv(eventDataFrame_wes, './results/mut_signatures/input_signal.csv', quote = F, row.names = F)

