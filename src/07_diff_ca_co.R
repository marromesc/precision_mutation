library(maftools)

setwd('/mnt/albyn/maria/precision_mutation')


# Load sample metadata ----------------------------------------------------

#wes
samples_wes <- as.data.frame(fread('./data/WES/DCIS_Precision_WES_All_Samples.txt'))
samples_wes <- samples_wes[samples_wes$qc_normal == 'pass' & samples_wes$qc_pdcis == 'pass',]
samples_wes <- samples_wes[samples_wes$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA)
                           & samples_wes$surgery_final == 'BCS',]
patients_wes <- unique(samples_wes$patient_id); length(patients_wes)


#nki panel
samples_nki_panel <- as.data.frame(fread('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_Panel_NKI_Samples.txt'))
samples_nki_panel <- samples_nki_panel[samples_nki_panel$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA)
                                       & samples_nki_panel$surgery_final == 'BCS',]
patients_nki_panel <- unique(samples_nki_panel$patient_id); length(patients_nki_panel)
patients_nki_panel <- patients_nki_panel[!(patients_nki_panel %in% patients_wes)]
samples_nki_panel <- samples_nki_panel[samples_nki_panel$patient_id %in% patients_nki_panel,]

#kcl panel
samples_kcl_panel <- as.data.frame(fread('./data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_Panel_Sloane_Samples.txt'))
samples_kcl_panel <- samples_kcl_panel[samples_kcl_panel$first_subseq_event %in% c('death', 'NA', 'ipsilateral IBC', NA)
                                       & samples_nki_panel$surgery_final == 'BCS',]
patients_kcl_panel <- unique(samples_kcl_panel$patient_id); length(patients_kcl_panel)
patients_kcl_panel <- patients_kcl_panel[!(patients_kcl_panel %in% patients_wes)]
samples_kcl_panel <- samples_kcl_panel[samples_kcl_panel$patient_id %in% patients_kcl_panel,]

patients_all <- unique(c(patients_wes, patients_nki_panel, patients_kcl_panel))
samples_all <- rbind(samples_wes[,colnames(samples_nki_panel)], samples_nki_panel, samples_kcl_panel[,colnames(samples_nki_panel)])


# Load filtered mutations -------------------------------------------------

#wes
eventDataFrame_indel <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Pindel_Filtered.rds')
eventDataFrame_mutect <- readRDS('./data/WES/DCIS_Precision_CaCo_WES_Mutect_Filtered.rds')

eventDataFrame_mutect <- eventDataFrame_mutect[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')]
eventDataFrame_indel <- eventDataFrame_indel[,c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                                'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                                'er', 'her2', 'grade', 'radiotherapy', 'batch')]

eventDataFrame_all <- rbind(eventDataFrame_mutect, eventDataFrame_indel)
eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_wes,]

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_wes <- eventDataFrame_all

#kcl panel
eventDataFrame_mutect <- readRDS('./data/Panel/DCIS_Precision_Panel_KCL/DCIS_Precision_CaCo_Panel_Sloane_Mutect_Filtered.rds')

eventDataFrame_all <- eventDataFrame_mutect[,c('ExonicFunc.refGene', 'Chr', 'Start', 'End', 'Ref', 'Alt',
                                               'AF_PDCIS', 'Gene.refGene', 'patient_id', 'first_subseq_event',
                                               'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_kcl_panel & !(eventDataFrame_all$patient_id %in% patients_wes),]
eventDataFrame_all$batch <- 'SLO_Panel'

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_kcl_panel,]

eventDataFrame_kcl_panel <- eventDataFrame_all

#nki panel
eventDataFrame_mutect <- readRDS('./data/Panel/DCIS_Precision_Panel_NKI/DCIS_Precision_CaCo_Panel_NKI_Mutect_Filtered.rds')

eventDataFrame_all <- eventDataFrame_mutect[,c('Consequence', 'CHROM', 'POS', 'POS', 'REF', 'ALT',
                                               'vaf_DCIS_DNA', 'Gene.refGene', 'patient_id', 'first_subseq_event',
                                               'er', 'her2', 'grade', 'radiotherapy')]

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_nki_panel & !(eventDataFrame_all$patient_id %in% patients_wes),]
eventDataFrame_all$batch <- 'NKI_Panel'

colnames(eventDataFrame_all) <- c('exonicfunc.knowngene', 'chr', 'start', 'end', 'ref_allele', 'alt_allele',
                                  'tumor_f', 'gene.knowngene', 'patient_id', 'first_subseq_event',
                                  'er', 'her2', 'grade', 'radiotherapy', 'batch')

eventDataFrame_all$case_control <- ifelse(eventDataFrame_all$first_subseq_event %in% c('ipsilateral IBC'), 'case',
                                          ifelse(eventDataFrame_all$first_subseq_event %in% c('NA', 'death'), 'control', 'dcis_case'))

eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients_nki_panel,]

eventDataFrame_nki_panel <- eventDataFrame_all

#rbind
eventDataFrame <- rbind(eventDataFrame_kcl_panel, eventDataFrame_nki_panel, eventDataFrame_wes)
eventDataFrame_panel <- rbind(eventDataFrame_nki_panel, eventDataFrame_kcl_panel)

# Differential analysis ---------------------------------------------------

for (eventDataFrame_all in list(eventDataFrame, eventDataFrame_wes, eventDataFrame_panel)){
  eventDataFrame_control <- eventDataFrame_all[eventDataFrame_all$case_control == 'control',]
  eventDataFrame_case <- eventDataFrame_all[eventDataFrame_all$case_control == 'case',]
  
  # sample size
  m2.sampleSize <- length(unique(eventDataFrame_control$patient_id))
  m1.sampleSize <- length(unique(eventDataFrame_case$patient_id))
  
  # get gene summary
  ct_case <- table(eventDataFrame_case$gene.knowngene)
  ct_control <- table(eventDataFrame_control$gene.knowngene)
  
  # remove mutations affecting less than 5 samples
  m1.genes = names(ct_case[ct_case > 3])
  m2.genes = names(ct_control[ct_control > 3])
  uniqueGenes = unique(c(m1.genes, m2.genes))
  
  # make frequency table
  m1.gs.comGenes = as.data.frame(ct_case[names(ct_case) %in% uniqueGenes])
  m2.gs.comGenes = as.data.frame(ct_control[names(ct_control) %in% uniqueGenes])
  
  # merge groups
  m.gs.meged = merge(m1.gs.comGenes, m2.gs.comGenes,
                     by = 'Var1', all = TRUE)
  
  #Set missing genes to zero
  m.gs.meged[is.na(m.gs.meged)] = 0
  m.gs.meged = as.data.frame(m.gs.meged)
  
  # get mut count for panel genes
  GenesPanel <- readRDS('./data/GenesPanel.RDS')
  m.gs.meged = m.gs.meged[m.gs.meged$Var1 %in% GenesPanel,]
  
  # fisher table
  fisherTable <- lapply(seq_len(nrow(m.gs.meged)), function(i){
    gene = m.gs.meged[i, 1]
    m1Mut = m.gs.meged[i,2]
    m2Mut = m.gs.meged[i,3]
    #print(i)
    
    ft_mat = matrix(c(m1Mut, m1.sampleSize-m1Mut, m2Mut, m2.sampleSize-m2Mut),
                    byrow = TRUE, nrow = 2)
    
    xf = fisher.test(ft_mat, conf.int = TRUE, conf.level = 0.95)
    
    pval = xf$p.value
    or = xf$estimate
    ci.up = xf$conf.int[2]
    ci.low = xf$conf.int[1]
    tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut , m2Mut,
                                  pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
    tdat
  })
  
  fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, fill = TRUE)
  fisherTable = fisherTable[order(pval)]
  
  # fdr
  fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]
  colnames(fisherTable)[2:3] = c('case', 'control')
  
  write.table(as.data.frame(fisherTable), paste0('./results/diff_mut_', paste(unique(eventDataFrame_all$batch), collapse = '_'), '.csv'), sep = ',', row.names = F)
}


# Differential analysis RT ---------------------------------------------------

for (rt in c(0, 1)){
  message(rt)
  eventDataFrame_all <- eventDataFrame
  patients <- samples_all$patient_id[samples_all$radiotherapy == rt]
  
  eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients,]
  
  eventDataFrame_control <- eventDataFrame_all[eventDataFrame_all$case_control == 'control',]
  eventDataFrame_case <- eventDataFrame_all[eventDataFrame_all$case_control == 'case',]
  
  # sample size
  m2.sampleSize <- length(unique(eventDataFrame_control$patient_id))
  m1.sampleSize <- length(unique(eventDataFrame_case$patient_id))
  
  # get gene summary
  ct_case <- table(eventDataFrame_case$gene.knowngene)
  ct_control <- table(eventDataFrame_control$gene.knowngene)
  
  # remove mutations affecting less than 5 samples
  m1.genes = names(ct_case[ct_case > 3])
  m2.genes = names(ct_control[ct_control > 3])
  uniqueGenes = unique(c(m1.genes, m2.genes))
  
  # make frequency table
  m1.gs.comGenes = as.data.frame(ct_case[names(ct_case) %in% uniqueGenes])
  m2.gs.comGenes = as.data.frame(ct_control[names(ct_control) %in% uniqueGenes])
  
  # merge groups
  m.gs.meged = merge(m1.gs.comGenes, m2.gs.comGenes,
                     by = 'Var1', all = TRUE)
  
  #Set missing genes to zero
  m.gs.meged[is.na(m.gs.meged)] = 0
  m.gs.meged = as.data.frame(m.gs.meged)
  
  # get mut count for panel genes
  GenesPanel <- readRDS('./data/GenesPanel.RDS')
  m.gs.meged = m.gs.meged[m.gs.meged$Var1 %in% GenesPanel,]
  
  # fisher table
  fisherTable <- lapply(seq_len(nrow(m.gs.meged)), function(i){
    gene = m.gs.meged[i, 1]
    m1Mut = m.gs.meged[i,2]
    m2Mut = m.gs.meged[i,3]
    #print(i)
    
    ft_mat = matrix(c(m1Mut, m1.sampleSize-m1Mut, m2Mut, m2.sampleSize-m2Mut),
                    byrow = TRUE, nrow = 2)
    
    xf = fisher.test(ft_mat, conf.int = TRUE, conf.level = 0.95)
    
    pval = xf$p.value
    or = xf$estimate
    ci.up = xf$conf.int[2]
    ci.low = xf$conf.int[1]
    tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut , m2Mut,
                                  pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
    tdat
  })
  
  fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, fill = TRUE)
  fisherTable = fisherTable[order(pval)]
  
  # fdr
  fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]
  colnames(fisherTable)[2:3] = c('case', 'control')
  
  write.table(as.data.frame(fisherTable), paste0('./results/diff_mut_RT_', rt, '.csv'), sep = ',', row.names = F)
}


# Differential analysis Her2 ---------------------------------------------------

eventDataFrame_all <- eventDataFrame
for (her2 in c(0, 1)){
  eventDataFrame_all <- eventDataFrame
  patients <- samples_all$patient_id[samples_all$her2 == her2]
  
  eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients,]
  
  eventDataFrame_control <- eventDataFrame_all[eventDataFrame_all$case_control == 'control',]
  eventDataFrame_case <- eventDataFrame_all[eventDataFrame_all$case_control == 'case',]
  
  # sample size
  m2.sampleSize <- length(unique(eventDataFrame_control$patient_id))
  m1.sampleSize <- length(unique(eventDataFrame_case$patient_id))
  
  # get gene summary
  ct_case <- table(eventDataFrame_case$gene.knowngene)
  ct_control <- table(eventDataFrame_control$gene.knowngene)
  
  # remove mutations affecting less than 5 samples
  m1.genes = names(ct_case[ct_case > 3])
  m2.genes = names(ct_control[ct_control > 3])
  uniqueGenes = unique(c(m1.genes, m2.genes))
  
  # make frequency table
  m1.gs.comGenes = as.data.frame(ct_case[names(ct_case) %in% uniqueGenes])
  m2.gs.comGenes = as.data.frame(ct_control[names(ct_control) %in% uniqueGenes])
  
  # merge groups
  m.gs.meged = merge(m1.gs.comGenes, m2.gs.comGenes,
                     by = 'Var1', all = TRUE)
  
  #Set missing genes to zero
  m.gs.meged[is.na(m.gs.meged)] = 0
  m.gs.meged = as.data.frame(m.gs.meged)
  
  # get mut count for panel genes
  GenesPanel <- readRDS('./data/GenesPanel.RDS')
  m.gs.meged = m.gs.meged[m.gs.meged$Var1 %in% GenesPanel,]
  
  # fisher table
  fisherTable <- lapply(seq_len(nrow(m.gs.meged)), function(i){
    gene = m.gs.meged[i, 1]
    m1Mut = m.gs.meged[i,2]
    m2Mut = m.gs.meged[i,3]
    #print(i)
    
    ft_mat = matrix(c(m1Mut, m1.sampleSize-m1Mut, m2Mut, m2.sampleSize-m2Mut),
                    byrow = TRUE, nrow = 2)
    
    xf = fisher.test(ft_mat, conf.int = TRUE, conf.level = 0.95)
    
    pval = xf$p.value
    or = xf$estimate
    ci.up = xf$conf.int[2]
    ci.low = xf$conf.int[1]
    tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut , m2Mut,
                                  pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
    tdat
  })
  
  fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, fill = TRUE)
  fisherTable = fisherTable[order(pval)]
  
  # fdr
  fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]
  colnames(fisherTable)[2:3] = c('case', 'control')
  
  write.table(as.data.frame(fisherTable), paste0('./results/diff_mut_HER2_', her2, '.csv'), sep = ',', row.names = F)
}

# Differential analysis ER ---------------------------------------------------

eventDataFrame_all <- eventDataFrame
for (er in c(0, 1)){
  eventDataFrame_all <- eventDataFrame
  patients <- samples_all$patient_id[samples_all$er == er]
  
  eventDataFrame_all <- eventDataFrame_all[eventDataFrame_all$patient_id %in% patients,]
  
  eventDataFrame_control <- eventDataFrame_all[eventDataFrame_all$case_control == 'control',]
  eventDataFrame_case <- eventDataFrame_all[eventDataFrame_all$case_control == 'case',]
  
  # sample size
  m2.sampleSize <- length(unique(eventDataFrame_control$patient_id))
  m1.sampleSize <- length(unique(eventDataFrame_case$patient_id))
  
  # get gene summary
  ct_case <- table(eventDataFrame_case$gene.knowngene)
  ct_control <- table(eventDataFrame_control$gene.knowngene)
  
  # remove mutations affecting less than 5 samples
  m1.genes = names(ct_case[ct_case > 3])
  m2.genes = names(ct_control[ct_control > 3])
  uniqueGenes = unique(c(m1.genes, m2.genes))
  
  # make frequency table
  m1.gs.comGenes = as.data.frame(ct_case[names(ct_case) %in% uniqueGenes])
  m2.gs.comGenes = as.data.frame(ct_control[names(ct_control) %in% uniqueGenes])
  
  # merge groups
  m.gs.meged = merge(m1.gs.comGenes, m2.gs.comGenes,
                     by = 'Var1', all = TRUE)
  
  #Set missing genes to zero
  m.gs.meged[is.na(m.gs.meged)] = 0
  m.gs.meged = as.data.frame(m.gs.meged)
  
  # get mut count for panel genes
  GenesPanel <- readRDS('./data/GenesPanel.RDS')
  m.gs.meged = m.gs.meged[m.gs.meged$Var1 %in% GenesPanel,]
  
  # fisher table
  fisherTable <- lapply(seq_len(nrow(m.gs.meged)), function(i){
    gene = m.gs.meged[i, 1]
    m1Mut = m.gs.meged[i,2]
    m2Mut = m.gs.meged[i,3]
    #print(i)
    
    ft_mat = matrix(c(m1Mut, m1.sampleSize-m1Mut, m2Mut, m2.sampleSize-m2Mut),
                    byrow = TRUE, nrow = 2)
    
    xf = fisher.test(ft_mat, conf.int = TRUE, conf.level = 0.95)
    
    pval = xf$p.value
    or = xf$estimate
    ci.up = xf$conf.int[2]
    ci.low = xf$conf.int[1]
    tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut , m2Mut,
                                  pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
    tdat
  })
  
  fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, fill = TRUE)
  fisherTable = fisherTable[order(pval)]
  
  # fdr
  fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]
  colnames(fisherTable)[2:3] = c('case', 'control')
  
  write.table(as.data.frame(fisherTable), paste0('./results/diff_mut_ER_', er, '.csv'), sep = ',', row.names = F)
}

