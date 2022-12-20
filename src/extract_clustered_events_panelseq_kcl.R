#read annotated for the annotations

# read the non filter to keep to take back the hostpot mutations which were 
# filtered out bc clustered events or other kind of potential artifact

filtered_vcf_files <- list.files('./data/DCIS_Precision_Panel_KCL', pattern = '*filtered.vcf', full.names = T)

samples <- gsub("\\_.*","",gsub(".*/","",filtered_vcf_files))
sloane <- as.data.frame(sloane[!(sloane$Tissue %in% c('Normal', 'DCIS recurrence', 'Inv recurrence', 'Inv recurrence Tubular', 'Inv recurrence NST', 'Inv recurrence 2')) 
                               & !is.na(sloane$'ng SureSelect'),])
rownames(sloane) <- sloane$'Sloane ID'
sloane <- sloane[samples,]

sloane <- merge(sloane, openclinica, by.x = 'Precision ID', by.y = 'patient_id')
Ind <- which(sloane$first_subseq_event %in% c('NA', 'ipsilateral DCIS', 'ipsilateral IBC'))

data_panel_kcl <- sloane[Ind,]

filtered_vcf <- lapply(filtered_vcf_files[Ind], readr::read_tsv, comment = "#", col_names = F)
names(filtered_vcf) <- sloane[Ind, 'DNA-Short']

toremove <- c()
for (i in 1:length(filtered_vcf)){
  message(i)
  if (nrow(filtered_vcf[[i]]) == 0){
    toremove <- c(toremove,i)
  } else {
    filtered_vcf[[i]] <- as.data.frame(filtered_vcf[[i]])
    filtered_vcf[[i]][,3] <- names(filtered_vcf)[i]
    colnames(filtered_vcf[[i]]) <- c('Chromosome', 'Position','ID','Ref','Alt','QUAL','Filter','Info','Format','Normal', 'PDCIS')
  }
}

filtered_vcf <- filtered_vcf[-toremove]
filtered_vcf <- do.call(rbind, filtered_vcf)

filtered_vcf$Info <- gsub('CONTQ=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('DP=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('ECNT=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('GERMQ=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('MBQ=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('MFRL=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('MMQ=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('MPOS=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('NALOD=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('NLOD=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('POPAF=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('SAAF=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('SAPP=','',filtered_vcf$Info)
filtered_vcf$Info <- gsub('TLOD=','',filtered_vcf$Info)
filtered_vcf <- filtered_vcf %>% separate(Info, sep=';', into = c('CONTQ', 'DP', 'ECNT', 'GERMQ', 'MBQ', 'MFRL', 'MMQ', 'MPOS', 'NALOD', 'NLOD', 'POPAF', 'SAAF', 'SAPP', 'TLOD'))
filtered_vcf <- filtered_vcf %>% separate(Normal, sep=':', into = c('GT_NOR', 'AD_NOR', 'AF_NOR', 'DP_NOR', 'F1R2_NOR', 'F2R2_NOR'))
filtered_vcf <- filtered_vcf %>% separate(PDCIS, sep=':', into = c('GT_PDCIS', 'AD_PDCIS', 'AF_PDCIS', 'DP_PDCIS', 'F1R2_PDCIS', 'F2R2_PDCIS'))
filtered_vcf <- filtered_vcf[,-8]

clustered_events <- filtered_vcf[grep('clustered_events', filtered_vcf$Filter),]

gene_annot <- readRDS('/mnt/albyn/maria/copynumber/precision_copy_number/gene_annot.RDS')
gene_annot <- gene_annot[gene_annot$chromosome_name %in% c(1:22, 'X') & gene_annot$hgnc_symbol %in% GenePanel,]

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
gr <- GRanges( seqnames = gene_annot$chromosome_name, ranges = IRanges(gene_annot$start_position, end = gene_annot$end_position), gene = gene_annot$hgnc_symbol)

seqlevelsStyle(gr) = "UCSC"  # necessary
cur19 = liftOver(gr, ch)
cur19 <- as.data.frame(cur19)

clustered_events$gene <- NA
for(i in 1:nrow(clustered_events)){
  res <- cur19[cur19$seqnames == clustered_events$Chromosome[i] & cur19$start <= clustered_events$Position[i] & cur19$end > clustered_events$Position[i],'gene']
  
  if (length(res) != 0){
    clustered_events$gene[i] <- res
  }
}

clustered_events_only <- clustered_events[clustered_events$Filter == 'clustered_events' & clustered_events$AF_PDCIS > 0.05,]

write.table(clustered_events, './results/clustered_events_Sloane_Target.txt', row.names = F, quote = F)
write.table(as.data.frame(table(clustered_events_only$gene)), './results/clustered_events_Sloane_Target_count_only_clustered_events.txt', row.names = F, quote = F)
write.table(as.data.frame(table(clustered_events$gene)), './results/clustered_events_Sloane_Target_count.txt', row.names = F, quote = F)

