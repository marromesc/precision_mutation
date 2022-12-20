GenesPanel <- readRDS('./data/GenesPanel.RDS')


# KCL ---------------------------------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

system('mkdir ./results/match_panel_wes')

panel_samples <- data.table::fread('./data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Samples.txt')
panel_samples <- panel_samples[panel_samples$first_subseq_event %in% c(NA, 'ipsilateral DCIS', 'ipsilateral IBC', 'death'),]

panel_filt_mut <- readRDS('./data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect_Filtered.rds'); dim(panel_filt_mut)

panel_unfilt_mut <- as.data.frame(data.table::fread('./data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect.txt'))
panel_unfilt_mut <- panel_unfilt_mut %>% dplyr::filter(Func.refGene %in% c('exonic', 'splicing', 'exonic;splicing'))
panel_unfilt_mut <- panel_unfilt_mut[!panel_unfilt_mut$ExonicFunc.refGene %in% c("synonymous SNV", 'unknown'),]
panel_unfilt_mut$ExonicFunc.refGene[panel_unfilt_mut$ExonicFunc.refGene %in% c(".")] <- "splicing"
panel_unfilt_mut$ExonicFunc.refGene[panel_unfilt_mut$ExonicFunc.refGene %in% c("nonsynonymous SNV")] <- "missense"
panel_unfilt_mut$ExonicFunc.refGene[panel_unfilt_mut$ExonicFunc.refGene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"
panel_unfilt_mut$ExonicFunc.refGene[panel_unfilt_mut$ExonicFunc.refGene %in% c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution")] <- "inframe_indel"
panel_unfilt_mut$ExonicFunc.refGene[panel_unfilt_mut$ExonicFunc.refGene %in% c("frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift"

wes_samples <- data.table::fread('./data/WES/DCIS_Precision_WES_All_Samples.txt'); dim(wes_samples)
wes_samples <- wes_samples[wes_samples$first_subseq_event %in% c(NA, 'ipsilateral DCIS', 'ipsilateral IBC', 'death'),]; dim(wes_samples)
wes_samples <- wes_samples[wes_samples$qc_pdcis == 'pass' & wes_samples$qc_normal == 'pass',]; dim(wes_samples)

eventDataFrame_mutect <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds')
eventDataFrame_indel <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds')

eventDataFrame_mutect <- eventDataFrame_mutect[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                  'tumor_f', 'gene.knowngene', 'batch', 'patient_id', 'first_subseq_event',
                                                  'er', 'her2', 'grade', 'radiotherapy')]
eventDataFrame_indel <- eventDataFrame_indel[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                'tumor_f', 'gene.knowngene', 'batch', 'patient_id', 'first_subseq_event',
                                                'er', 'her2', 'grade', 'radiotherapy')]

wes_filt_mut <- as.data.frame(rbind(eventDataFrame_mutect, eventDataFrame_indel))

wes_unfilt_mut <- as.data.frame(readRDS('./data/WES/DCIS_Precision_CaCo_WES_Mutect.rds')); dim(wes_unfilt_mut)
wes_unfilt_mut <- wes_unfilt_mut %>% dplyr::filter(func.knowngene %in% c('exonic', 'splicing', 'exonic;splicing')); dim(wes_unfilt_mut)
wes_unfilt_mut <- wes_unfilt_mut[!wes_unfilt_mut$exonicfunc.knowngene %in% c("synonymous SNV", 'unknown'),]; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c(".")] <- "splicing"; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c("nonsynonymous SNV")] <- "missense"; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution")] <- "inframe_indel"; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c("frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift"; dim(wes_unfilt_mut)

matched_samples <- panel_samples$patient_id[panel_samples$patient_id %in% wes_samples$patient_id]

wes_unfilt_mut$mut_code <- paste0(wes_unfilt_mut$gene.knowngene, '-', wes_unfilt_mut$chr, ':', as.integer(wes_unfilt_mut$start), '-', wes_unfilt_mut$ref_allele, '>', wes_unfilt_mut$alt_allele, '-', wes_unfilt_mut$exonicfunc.knowngene, '-', wes_unfilt_mut$precision_patient_id) 
wes_filt_mut$mut_code <- paste0(wes_filt_mut$gene.knowngene, '-', wes_filt_mut$chr, ':', as.integer(wes_filt_mut$start), '-', wes_filt_mut$ref_allele, '>', wes_filt_mut$alt_allele, '-', wes_filt_mut$exonicfunc.knowngene, '-', wes_filt_mut$precision_patient_id) 
panel_unfilt_mut$mut_code <- paste0(panel_unfilt_mut$Gene.refGene, '-', panel_unfilt_mut$Chr, ':', as.integer(panel_unfilt_mut$Start), '-', panel_unfilt_mut$Ref, '>', panel_unfilt_mut$Alt, '-', panel_unfilt_mut$ExonicFunc.refGene, '-', panel_unfilt_mut$patient_id) 
panel_filt_mut$mut_code <- paste0(panel_filt_mut$Gene.refGene, '-', panel_filt_mut$Chr, ':', as.integer(panel_filt_mut$Start), '-', panel_filt_mut$Ref, '>', panel_filt_mut$Alt, '-', panel_filt_mut$ExonicFunc.refGen, '-', panel_filt_mut$patient_id) 

discrepant <- data.frame(mut_cod = unique(c(wes_unfilt_mut$mut_code[wes_unfilt_mut$gene.knowngene %in% GenesPanel & wes_unfilt_mut$precision_patient_id %in% matched_samples], wes_filt_mut$mut_code[wes_filt_mut$gene.knowngene %in% GenesPanel & wes_filt_mut$precision_patient_id %in% matched_samples], panel_unfilt_mut$mut_code[panel_unfilt_mut$Gene.refGene %in% GenesPanel & panel_unfilt_mut$patient_id %in% matched_samples], panel_filt_mut$mut_code[panel_filt_mut$Gene.refGene %in% GenesPanel & panel_filt_mut$patient_id %in% matched_samples])), gene = NA, sample = NA, chrom = NA, pos = NA, ref = NA, alt = NA, consequence = NA, vaf_panel = NA, vaf_wes = NA, filtered_panel = NA, filtered_wes = NA)
for (i in 1:nrow(discrepant)){
  message(i)
  discrepant$gene[i] <- strsplit(discrepant$mut_cod[i], '-')[[1]][1]
  discrepant$sample[i] <- strsplit(discrepant$mut_cod[i], '-')[[1]][5]
  discrepant$chrom[i] <- strsplit(strsplit(discrepant$mut_cod[i], '-')[[1]][2], ':')[[1]][1]
  discrepant$pos[i] <- strsplit(strsplit(discrepant$mut_cod[i], '-')[[1]][2], ':')[[1]][2]
  discrepant$ref[i] <- strsplit(strsplit(discrepant$mut_cod[i], '-')[[1]][3], '>')[[1]][1]
  discrepant$alt[i] <- strsplit(strsplit(discrepant$mut_cod[i], '-')[[1]][3], '>')[[1]][2]
  discrepant$consequence[i] <- strsplit(discrepant$mut_cod[i], '-')[[1]][4]
  
  Ind_unfilt_wes <- which(wes_unfilt_mut$mut_code == discrepant$mut_cod[i])
  Ind_filt_wes <- which(wes_filt_mut$mut_code == discrepant$mut_cod[i])
  Ind_unfilt_panel <- which(panel_unfilt_mut$mut_code == discrepant$mut_cod[i])
  Ind_filt_panel <- which(panel_filt_mut$mut_code == discrepant$mut_cod[i])
  
  if ( length(Ind_unfilt_wes)>0 | length(Ind_filt_wes)>0 ){
    discrepant$vaf_wes[i] <- wes_unfilt_mut$tumor_f[Ind_unfilt_wes]
    if (is.na(discrepant$vaf_wes[i])){
      discrepant$vaf_wes[i] <- wes_filt_mut$tumor_f[Ind_filt_wes]
    }
  }
  
  if(length(Ind_filt_wes)>0){
    discrepant$filtered_wes[i] <- 'Pass'
  } else {
    discrepant$filtered_wes[i] <- 'Fail'
  }
  
  if ( length(Ind_unfilt_panel)>0 | length(Ind_filt_panel)>0 ){
    discrepant$vaf_panel[i] <- panel_unfilt_mut$AF_PDCIS[Ind_unfilt_panel]
    if (is.na(discrepant$vaf_panel[i])){
      discrepant$vaf_panel[i] <- panel_filt_mut$AF_PDCIS[Ind_filt_panel]
    }
  }
  
  if(length(Ind_filt_panel)>0){
    discrepant$filtered_panel[i] <- 'Pass'
  } else {
    discrepant$filtered_panel[i] <- 'Fail'
  }
}

write.csv(discrepant, './results/match_panel_wes/discrepant_muts.csv')

  for (pat in matched_samples){
  message(pat)
  panel_sample <- panel_filt_mut[panel_filt_mut$patient_id == pat,]
  wes_sample <- wes_filt_mut[wes_filt_mut$patient_id == pat,]
  panel_sample_unf <- panel_unfilt_mut[panel_unfilt_mut$patient_id == pat,]
  wes_sample_unf <- wes_unfilt_mut[wes_unfilt_mut$precision_patient_id == pat,]
  
  # prepare matrix
  GenesPanel <- readRDS('./data/GenesPanel.RDS')
  geneMatrix <- matrix("",nrow=4,ncol=length(GenesPanel))
  rownames(geneMatrix) <- c(paste0(pat, '_panel_unfiltered'), paste0(pat, '_panel'), paste0(pat, '_wes_unfiltered'), paste0(pat, '_wes'))
  colnames(geneMatrix) <- GenesPanel
 
  for (j in 1:ncol(geneMatrix)) {
    Ind_panel <- which(panel_sample$Gene.refGene == colnames(geneMatrix)[j])
    Ind_panel_unf <- which(panel_sample_unf$Gene.refGene == colnames(geneMatrix)[j])
    Ind_wes <- which(wes_sample$gene.knowngene == colnames(geneMatrix)[j])
    Ind_wes_unf <- which(wes_sample_unf$gene.knowngene == colnames(geneMatrix)[j])
    
    if (length(Ind_panel)==1){
      geneMatrix[paste0(pat, '_panel'),j] <- panel_sample$ExonicFunc.refGene[Ind_panel]
    } else if (length(Ind_panel)>1){
      geneMatrix[paste0(pat, '_panel'),j] <- "multi_hit"
    }
    
    if (length(Ind_wes_unf)==1){
      geneMatrix[paste0(pat, '_wes_unfiltered'),j] <- wes_sample_unf$exonicfunc.knowngene[Ind_wes_unf]
    } else if (length(Ind_wes_unf)>1){
      geneMatrix[paste0(pat, '_wes_unfiltered'),j] <- "multi_hit"
    }

    if (length(Ind_wes)==1){
      geneMatrix[paste0(pat, '_wes'),j] <- wes_sample$exonicfunc.knowngene[Ind_wes]
    } else if (length(Ind_wes)>1){
      geneMatrix[paste0(pat, '_wes'),j] <- "multi_hit"
    }
    
    if (length(Ind_panel_unf)==1){
      geneMatrix[paste0(pat, '_panel_unfiltered'),j] <- panel_sample_unf$ExonicFunc.refGene[Ind_panel_unf]
    } else if (length(Ind_panel_unf)>1){
      geneMatrix[paste0(pat, '_panel_unfiltered'),j] <- "multi_hit"
    }
  }
  
  pdf(paste0("./results/match_panel_wes/oncoPrint_", pat,".pdf"), height = 5)
  print(oncoPrint(t(geneMatrix),
                  col=col,
                  alter_fun=alter_fun,
                  show_column_names=TRUE,
                  remove_empty_columns = FALSE,
                  remove_empty_rows = FALSE,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 5),
                  show_pct=TRUE,
                  column_order=rownames(geneMatrix),
                  row_order=genes_sorted,
                  #top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot()),
                  #right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE)),
                  heatmap_legend_param = list(title = "Alterations"),
                  #bottom_annotation = ha,
                  alter_fun_is_vectorized = FALSE,
                  pct_gp = gpar(fontsize = 6),
                  #column_title = paste0(unique(eventDataFrame$case_control), ' in Panel-Seq (n=', length(patients), ')')
  ))
  dev.off()
}



# NKI ---------------------------------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

system('mkdir ./results/match_panel_wes')

panel_samples <- data.table::fread('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples.txt')
panel_samples <- panel_samples[panel_samples$first_subseq_event %in% c(NA, 'ipsilateral DCIS', 'ipsilateral IBC', 'death'),]

panel_filt_mut <- readRDS('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.rds'); dim(panel_filt_mut)

panel_unfilt_mut <- as.data.frame(readxl::read_xlsx('./data/Panel/DCIS_Precision_Panel_NKI/Compiled_mutations_nki.xlsx'))
panel_unfilt_mut$Consequence[panel_unfilt_mut$Consequence %in% c("stop_gained","stoploss","startgain","startloss")] <- "nonsense"; dim(panel_unfilt_mut)
panel_unfilt_mut$Consequence[panel_unfilt_mut$Consequence %in% c('splice_donor_variant&intron_variant',"missense_variant", "missense_variant&splice_region_variant")] <- "missense"; dim(panel_unfilt_mut)
panel_unfilt_mut$Consequence[panel_unfilt_mut$Consequence %in% c('frameshift_variant&splice_region_variant',"frameshift_variant", "frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift"; dim(panel_unfilt_mut)
panel_unfilt_mut$Consequence[panel_unfilt_mut$Consequence %in% c(".", 'splice_donor_variant', 'splice_acceptor_variant')] <- "splicing"; dim(panel_unfilt_mut)

wes_samples <- data.table::fread('./data/WES/DCIS_Precision_WES_All_Samples.txt'); dim(wes_samples)
wes_samples <- wes_samples[wes_samples$first_subseq_event %in% c(NA, 'ipsilateral DCIS', 'ipsilateral IBC', 'death'),]; dim(wes_samples)
wes_samples <- wes_samples[wes_samples$qc_pdcis == 'pass' & wes_samples$qc_normal == 'pass',]; dim(wes_samples)

eventDataFrame_mutect <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds')
eventDataFrame_indel <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds')

eventDataFrame_mutect <- eventDataFrame_mutect[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                  'tumor_f', 'gene.knowngene', 'batch', 'patient_id', 'first_subseq_event',
                                                  'er', 'her2', 'grade', 'radiotherapy')]
eventDataFrame_indel <- eventDataFrame_indel[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                'tumor_f', 'gene.knowngene', 'batch', 'patient_id', 'first_subseq_event',
                                                'er', 'her2', 'grade', 'radiotherapy')]

wes_filt_mut <- as.data.frame(rbind(eventDataFrame_mutect, eventDataFrame_indel))

wes_unfilt_mut <- as.data.frame(readRDS('./data/WES/DCIS_Precision_CaCo_WES_Mutect.rds')); dim(wes_unfilt_mut)
wes_unfilt_mut <- wes_unfilt_mut %>% dplyr::filter(func.knowngene %in% c('exonic', 'splicing', 'exonic;splicing')); dim(wes_unfilt_mut)
wes_unfilt_mut <- wes_unfilt_mut[!wes_unfilt_mut$exonicfunc.knowngene %in% c("synonymous SNV", 'unknown'),]; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c(".")] <- "splicing"; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c("nonsynonymous SNV")] <- "missense"; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution")] <- "inframe_indel"; dim(wes_unfilt_mut)
wes_unfilt_mut$exonicfunc.knowngene[wes_unfilt_mut$exonicfunc.knowngene %in% c("frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift"; dim(wes_unfilt_mut)

matched_samples <- panel_samples$patient_id[panel_samples$patient_id %in% wes_samples$patient_id]

wes_unfilt_mut$mut_code <- paste0(wes_unfilt_mut$gene.knowngene, '-', wes_unfilt_mut$chr, ':', as.integer(wes_unfilt_mut$start), '-', wes_unfilt_mut$ref_allele, '>', wes_unfilt_mut$alt_allele, '-', wes_unfilt_mut$exonicfunc.knowngene, '-', wes_unfilt_mut$precision_patient_id) 
wes_filt_mut$mut_code <- paste0(wes_filt_mut$gene.knowngene, '-', wes_filt_mut$chr, ':', as.integer(wes_filt_mut$start), '-', wes_filt_mut$ref_allele, '>', wes_filt_mut$alt_allele, '-', wes_filt_mut$exonicfunc.knowngene, '-', wes_filt_mut$precision_patient_id) 
panel_unfilt_mut$mut_code <- paste0(panel_unfilt_mut$Gene.refGene, '-', panel_unfilt_mut$Chr, ':', as.integer(panel_unfilt_mut$Start), '-', panel_unfilt_mut$Ref, '>', panel_unfilt_mut$Alt, '-', panel_unfilt_mut$ExonicFunc.refGene, '-', panel_unfilt_mut$patient_id) 
panel_filt_mut$mut_code <- paste0(panel_filt_mut$Gene.refGene, '-', panel_filt_mut$Chr, ':', as.integer(panel_filt_mut$Start), '-', panel_filt_mut$Ref, '>', panel_filt_mut$Alt, '-', panel_filt_mut$ExonicFunc.refGen, '-', panel_filt_mut$patient_id) 

discrepant <- data.frame(mut_cod = unique(c(wes_unfilt_mut$mut_code[wes_unfilt_mut$gene.knowngene %in% GenesPanel & wes_unfilt_mut$precision_patient_id %in% matched_samples], wes_filt_mut$mut_code[wes_filt_mut$gene.knowngene %in% GenesPanel & wes_filt_mut$precision_patient_id %in% matched_samples], panel_unfilt_mut$mut_code[panel_unfilt_mut$Gene.refGene %in% GenesPanel & panel_unfilt_mut$patient_id %in% matched_samples], panel_filt_mut$mut_code[panel_filt_mut$Gene.refGene %in% GenesPanel & panel_filt_mut$patient_id %in% matched_samples])), gene = NA, sample = NA, chrom = NA, pos = NA, ref = NA, alt = NA, consequence = NA, vaf_panel = NA, vaf_wes = NA, filtered_panel = NA, filtered_wes = NA)
for (i in 1:nrow(discrepant)){
  message(i)
  discrepant$gene[i] <- strsplit(discrepant$mut_cod[i], '-')[[1]][1]
  discrepant$sample[i] <- strsplit(discrepant$mut_cod[i], '-')[[1]][5]
  discrepant$chrom[i] <- strsplit(strsplit(discrepant$mut_cod[i], '-')[[1]][2], ':')[[1]][1]
  discrepant$pos[i] <- strsplit(strsplit(discrepant$mut_cod[i], '-')[[1]][2], ':')[[1]][2]
  discrepant$ref[i] <- strsplit(strsplit(discrepant$mut_cod[i], '-')[[1]][3], '>')[[1]][1]
  discrepant$alt[i] <- strsplit(strsplit(discrepant$mut_cod[i], '-')[[1]][3], '>')[[1]][2]
  discrepant$consequence[i] <- strsplit(discrepant$mut_cod[i], '-')[[1]][4]
  
  Ind_unfilt_wes <- which(wes_unfilt_mut$mut_code == discrepant$mut_cod[i])
  Ind_filt_wes <- which(wes_filt_mut$mut_code == discrepant$mut_cod[i])
  Ind_unfilt_panel <- which(panel_unfilt_mut$mut_code == discrepant$mut_cod[i])
  Ind_filt_panel <- which(panel_filt_mut$mut_code == discrepant$mut_cod[i])
  
  if ( length(Ind_unfilt_wes)>0 | length(Ind_filt_wes)>0 ){
    discrepant$vaf_wes[i] <- wes_unfilt_mut$tumor_f[Ind_unfilt_wes]
    if (is.na(discrepant$vaf_wes[i])){
      discrepant$vaf_wes[i] <- wes_filt_mut$tumor_f[Ind_filt_wes]
    }
  }
  
  if(length(Ind_filt_wes)>0){
    discrepant$filtered_wes[i] <- 'Pass'
  } else {
    discrepant$filtered_wes[i] <- 'Fail'
  }
  
  if ( length(Ind_unfilt_panel)>0 | length(Ind_filt_panel)>0 ){
    discrepant$vaf_panel[i] <- panel_unfilt_mut$AF_PDCIS[Ind_unfilt_panel]
    if (is.na(discrepant$vaf_panel[i])){
      discrepant$vaf_panel[i] <- panel_filt_mut$AF_PDCIS[Ind_filt_panel]
    }
  }
  
  if(length(Ind_filt_panel)>0){
    discrepant$filtered_panel[i] <- 'Pass'
  } else {
    discrepant$filtered_panel[i] <- 'Fail'
  }
}

write.csv(discrepant, './results/match_panel_wes/discrepant_muts.csv')

for (pat in matched_samples){
  message(pat)
  samp <- panel_samples$panel_id[panel_samples$patient_id == pat]
  panel_sample <- panel_filt_mut[panel_filt_mut$patient_id == pat,]
  wes_sample <- wes_filt_mut[wes_filt_mut$patient_id == pat,]
  panel_sample_unf <- panel_unfilt_mut[panel_unfilt_mut$sample_ID_DCIS == samp,]
  wes_sample_unf <- wes_unfilt_mut[wes_unfilt_mut$precision_patient_id == pat,]
  
  # prepare matrix
  GenesPanel <- readRDS('./data/GenesPanel.RDS')
  geneMatrix <- matrix("",nrow=4,ncol=length(GenesPanel))
  rownames(geneMatrix) <- c(paste0(pat, '_panel_unfiltered'), paste0(pat, '_panel'), paste0(pat, '_wes_unfiltered'), paste0(pat, '_wes'))
  colnames(geneMatrix) <- GenesPanel
  
  for (j in 1:ncol(geneMatrix)) {
    Ind_panel <- which(panel_sample$Gene.refGene == colnames(geneMatrix)[j])
    Ind_panel_unf <- which(panel_sample_unf$SYMBOL == colnames(geneMatrix)[j])
    Ind_wes <- which(wes_sample$gene.knowngene == colnames(geneMatrix)[j])
    Ind_wes_unf <- which(wes_sample_unf$gene.knowngene == colnames(geneMatrix)[j])
    
    if (length(Ind_panel)==1){
      geneMatrix[paste0(pat, '_panel'),j] <- panel_sample$Consequence[Ind_panel]
    } else if (length(Ind_panel)>1){
      geneMatrix[paste0(pat, '_panel'),j] <- "multi_hit"
    }
    
    if (length(Ind_wes_unf)==1){
      geneMatrix[paste0(pat, '_wes_unfiltered'),j] <- wes_sample_unf$exonicfunc.knowngene[Ind_wes_unf]
    } else if (length(Ind_wes_unf)>1){
      geneMatrix[paste0(pat, '_wes_unfiltered'),j] <- "multi_hit"
    }
    
    if (length(Ind_wes)==1){
      geneMatrix[paste0(pat, '_wes'),j] <- wes_sample$exonicfunc.knowngene[Ind_wes]
    } else if (length(Ind_wes)>1){
      geneMatrix[paste0(pat, '_wes'),j] <- "multi_hit"
    }
    
    if (length(Ind_panel_unf)==1){
      geneMatrix[paste0(pat, '_panel_unfiltered'),j] <- panel_sample_unf$Consequence[Ind_panel_unf]
    } else if (length(Ind_panel_unf)>1){
      geneMatrix[paste0(pat, '_panel_unfiltered'),j] <- "multi_hit"
    }
  }
  
  pdf(paste0("./results/match_panel_wes/oncoPrint_", pat,".pdf"), height = 5)
  print(oncoPrint(t(geneMatrix),
                  col=col,
                  alter_fun=alter_fun,
                  show_column_names=TRUE,
                  remove_empty_columns = FALSE,
                  remove_empty_rows = FALSE,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 5),
                  show_pct=TRUE,
                  column_order=rownames(geneMatrix),
                  row_order=genes_sorted,
                  #top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot()),
                  #right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE)),
                  heatmap_legend_param = list(title = "Alterations"),
                  #bottom_annotation = ha,
                  alter_fun_is_vectorized = FALSE,
                  pct_gp = gpar(fontsize = 6),
                  #column_title = paste0(unique(eventDataFrame$case_control), ' in Panel-Seq (n=', length(patients), ')')
  ))
  dev.off()
}
