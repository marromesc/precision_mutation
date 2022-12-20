# function to read mutect files and save as a df
read_mutect = function(datapath, source, pattern = NULL){
        
        # combine all mutects 
        if (is.null(pattern)){
                ls.files = list.files(path = c(datapath), pattern = '.tsv', full.names = TRUE)
        } else {
                ls.files = do.call(rbind, lapply(pattern, function(x) list.files(path = c(datapath),
                                                                                 pattern = x, full.names = TRUE)))[,1]
        }

        
        datlist = list()
        for(i in 1:length(ls.files)){
                message(ls.files[i])
                if( file.info(ls.files[i])$size > 0 ){
                        dat.mat = read.delim2(file = ls.files[i]) %>% 
                                tbl_df() %>% 
                                mutate(tool = "mutect", source = source)
                        
                        datlist[[i]] = dat.mat
                } 
        }
        
        df_comb = do.call(rbind, datlist);dim(df_comb) # 2536    161
        # save the file
        #write.csv(df_comb, paste0("df_comb_", source, ".csv"))
        return(df_comb)
        
}

filterMaf = function(maf, genes = NULL, tsb = NULL, isTCGA = FALSE){
        
        if(all(c(is.null(tsb), is.null(genes)))){
                stop("Please provide sample names or genes or a query or ranges to subset by.")
        }
        
        #Synonymous variants
        maf.silent <- maf@maf.silent
        #Main data
        maf.dat <- maf@data
        #Annotations
        maf.anno <- data.table::copy(x = maf@clinical.data)
        nrows_nsyn = nrow(maf.dat)
        nrows_syn = nrow(maf.silent)
        
        #Select
        if(!is.null(tsb)){
                #message("-subsetting by tumor sample barcodes..")
                tsb = as.character(tsb)
                if(isTCGA){
                        tsb = substr(x = tsb, start = 1, stop = 12)
                }
                maf.dat = maf.dat[!Tumor_Sample_Barcode %in% tsb,]
                maf.silent = maf.silent[!Tumor_Sample_Barcode %in% tsb,]
                message("Removed ", (nrows_syn+nrows_nsyn) - (nrow(maf.dat) + nrow(maf.silent)), " variants from ", length(tsb), " samples")
        }
        
        if(!is.null(genes)){
                #message("-subsetting by genes..")
                genes = as.character(genes)
                maf.dat = maf.dat[!Hugo_Symbol %in% genes, ]
                maf.silent = maf.silent[!Hugo_Symbol %in% genes, ]
                message("Removed ", (nrows_syn+nrows_nsyn) - (nrow(maf.dat) + nrow(maf.silent)), " variants from ", length(genes), " genes")
        }
        
        maf.silent = droplevels.data.frame(maf.silent)
        maf.dat = droplevels.data.frame(maf.dat)
        maf.anno = droplevels.data.frame(maf.anno)
        
        MAF(nonSyn = maf.dat, syn = maf.silent, clinicalData = maf.anno, verbose = FALSE)
        
        # mafSummary = summarizeMaf(maf.dat, chatty = FALSE, anno = maf.anno)
        #
        # MAF(data = maf.dat, variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
        #         variant.classification.summary = mafSummary$variant.classification.summary, gene.summary = mafSummary$gene.summary,
        #         summary = mafSummary$summary, maf.silent = maf.silent, clinical.data = droplevels(mafSummary$sample.anno))
}
