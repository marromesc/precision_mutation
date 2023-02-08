#downloaded from cbioportal https://www.cbioportal.org/results/download?cancer_study_list=brca_tcga&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cgistic&case_set_id=brca_tcga_cnaseq&gene_list=TP53%252C%2520PIK3CA%252C%2520GATA3%252C%2520KMT2D%252C%2520NF1%252C%2520MAP3K1%252C%2520KMT2C%252C%2520ARID1A%252C%2520ATM%252C%2520ERBB3%252C%2520PIK3R1%252C%2520PTEN%252C%2520AKT1%252C%2520ARID1B%252C%2520CBFB%252C%2520BRCA1%252C%2520NOTCH1%252C%2520SF3B1%252C%2520TBX3%252C%2520BRCA2%252C%2520ERBB2%252C%2520RUNX1%252C%2520PDGFRA%252C%2520EGFR%252C%2520RB1%252C%2520NCOR1%252C%2520BRAF%252C%2520MYC%252C%2520IGF1R%252C%2520CDH1%252C%2520PBRM1%252C%2520ESR1%252C%2520FGFR2%252C%2520SMAD4%252C%2520STK11%252C%2520MAP2K4%252C%2520CHEK2%252C%2520FBXW7%252C%2520FGFR1%252C%2520CDKN2A%252C%2520CCND1%252C%2520BAP1%252C%2520CCND3%252C%2520MDM2%252C%2520CCNE1&geneset_list=%20&tab_index=tab_visualize&Action=Submit

dataset <- readr::read_tsv('alterations_across_samples_tcga.tsv')

geneMatrix <- as.matrix(dataset[,GenesPanel])
#colnames(geneMatrix) <- metabric$`Sample ID`
geneMatrix[geneMatrix == 'no alteration'] <- ''
geneMatrix[geneMatrix != ''] <- 'mutated'

Ind <- sample(1:nrow(geneMatrix), 980)
geneMatrix <- geneMatrix[Ind,]

col = c(mutated = "violet")
alter_fun <- list(
  mutated=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["mutated"], col=NA))
  }
)

pdf(paste0("oncoPrint_metabric.pdf"), height = 5)
print(oncoPrint(t(geneMatrix),
              col=col,
               alter_fun=alter_fun,
                show_column_names=F,
                remove_empty_columns = FALSE,
                remove_empty_rows = FALSE,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 5),
                show_pct=TRUE,
                column_order=rownames(geneMatrix),
                #row_order=genes_sorted,
                #top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot()),
                #right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE)),
                heatmap_legend_param = list(title = "Alterations"),
                #bottom_annotation = ha,
                alter_fun_is_vectorized = FALSE,
                pct_gp = gpar(fontsize = 6),
                #column_title = paste0(unique(eventDataFrame$case_control), ' in Panel-Seq (n=', length(patients), ')')
))
dev.off()
