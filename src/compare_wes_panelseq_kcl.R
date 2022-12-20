# Match WES and PanelSeq Sloane -------------------------------------------

bed <- read.delim('/mnt/albyn/maria/prj_precision/targetedseq/Sloane_Vandna/Sloane_Covered.bed', skip = 3, sep = "\t", header = FALSE, as.is = TRUE)
GenePanel <- unique(bed[,4])

wes_filtered <- read.delim("./DCIS_Precision_WES_Filtered.txt", sep = '\t', header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
wes_filtered <- wes_filtered[wes_filtered$Hugo_Symbol %in% GenePanel,]
panel_filtered <- read.delim("./DCIS_Precision_Panel_Sloane.txt", sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)

match_samples <- unique(wes_filtered$sample_label)
match_samples <- match_samples[match_samples %in% panel_filtered$ID]

wes_match <- wes_filtered[wes_filtered$sample_label %in% match_samples,]
wes_match$n_AF <- wes_match$n_alt_count / (wes_match$n_ref_count + wes_match$n_alt_count )
wes_match <- wes_match %>% dplyr::select('sample_label', 'Hugo_Symbol', 'Chromosome', 'Start_position',
                                         'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele2',
                                         't_alt_count', 't_ref_count', 'n_alt_count', 'n_ref_count', 'n_AF','VAF', 'exac_all',
                                         'esp6500siv2_all', 'ALL.sites.2015_08')
colnames(wes_match) <- c('id', 'gene', 'chr', 'pos',
                         'variant_classification', 'ref', 'alt',
                         't_alt_count', 't_ref_count', 'n_alt_count', 'n_ref_count', 'n_AF','VAF', 'exac_all',
                         'esp6500siv2_all', 'ALL.sites.2015_08')
panel_match <- panel_filtered[panel_filtered$ID %in% match_samples,]
panel_match$t_alt_count <- sapply(panel_match$AD_PDCIS, function(x) { x = strsplit(x, ',') ; as.numeric(x[[1]][2]) } )
panel_match$t_ref_count <- sapply(panel_match$AD_PDCIS, function(x) { x = strsplit(x, ',') ; as.numeric(x[[1]][1]) } )

panel_match$n_alt_count <- sapply(panel_match$AD_NOR, function(x) { x = strsplit(x, ',') ; as.numeric(x[[1]][2]) } )
panel_match$n_ref_count <- sapply(panel_match$AD_NOR, function(x) { x = strsplit(x, ',') ; as.numeric(x[[1]][1]) } )

panel_match <- panel_match %>% dplyr::select('ID', 'Gene.refGene', 'Chromosome', 'Position',
                                             'ExonicFunc.refGene', 'Ref', 'Alt',
                                             't_alt_count', 't_ref_count', 'n_alt_count', 'n_ref_count', 'AF_NOR', 'AF_PDCIS', 'ExAC_ALL',
                                             'esp6500siv2_all', 'ALL.sites.2015_08')
colnames(panel_match) <- c('id', 'gene', 'chr', 'pos',
                           'variant_classification', 'ref', 'alt',
                           't_alt_count', 't_ref_count', 'n_alt_count', 'n_ref_count', 'n_AF','VAF', 'exac_all',
                           'esp6500siv2_all', 'ALL.sites.2015_08')

wes_match$variant_classification[wes_match$variant_classification %in% c("Splice_Site")] <- "splicing"
wes_match$variant_classification[wes_match$variant_classification %in% c("Missense_Mutation")] <- "missense"
wes_match$variant_classification[wes_match$variant_classification %in% c("In_Frame_Del","In_Frame_Ins","nonframeshift_substitution")] <- "inframe_indel"
wes_match$variant_classification[wes_match$variant_classification %in% c("Frame_Shift_Del","Frame_Shift_Ins")] <- "frameshift"
wes_match$variant_classification[wes_match$variant_classification %in% c("Nonsense_Mutation","stoploss","startgain","startloss")] <- "nonsense"

panel_match$variant_classification[panel_match$variant_classification %in% c(".")] <- "splicing"
panel_match$variant_classification[panel_match$variant_classification %in% c("Missense_Mutation")] <- "missense"
panel_match$variant_classification[panel_match$variant_classification %in% c("nonframeshift_deletion","In_Frame_Ins","nonframeshift_substitution")] <- "inframe_indel"
panel_match$variant_classification[panel_match$variant_classification %in% c("frameshift_insertion","frameshift_deletion")] <- "frameshift"
panel_match$variant_classification[panel_match$variant_classification %in% c("Nonsense_Mutation","stopgain","stoploss","startloss")] <- "nonsense"

panel_match <- panel_match[panel_match$variant_classification %in%  c('splicing', 'missense', 'inframe_indel', 'frameshift', 'nonsense'),]

panel_match$chr <- gsub('chr','',panel_match$chr)

wes_match$mutation <- paste0(wes_match$id, wes_match$chr, ":", wes_match$pos, ':', wes_match$ref, '-', wes_match$alt, ':', wes_match$variant_classification)
panel_match$mutation <- paste0(panel_match$id, panel_match$chr, ":", panel_match$pos, ':', panel_match$ref, '-', panel_match$alt, ':', panel_match$variant_classification)

wes_match <- wes_match[order(wes_match$id),]
panel_match <- panel_match[order(panel_match$id),]

wes_match$mutation[wes_match$mutation %in% panel_match]

wes_match$mutation[wes_match$mutation %in% panel_match$mutation]

wes_match <- split(wes_match, wes_match$id)
panel_match <- split(panel_match, panel_match$id)

