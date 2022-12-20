vc_cols = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(vc_cols) = c(
    'Frame_Shift_Del',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Multi_Hit',
    'Frame_Shift_Ins',
    'In_Frame_Ins',
    'Splice_Site',
    'In_Frame_Del',
    'Gain',
    'Loss'
)

svc_cols = RColorBrewer::brewer.pal(n = 7, name = 'Paired')
names(svc_cols) = c(
    'Domestic_Ins',
    'Translocation',
    'Inversion',
    'Multi_Hit',
    'Foreign_Ins',
    'Duplication',
    'Deletion'
)

annotsv2MAF <- function(inputFile){
    # more columns than column names
    first_line <- read_lines(inputFile, n_max = 1)
    col_names <- unlist(strsplit(first_line, "\t"))
    tmp <- read.delim(inputFile, sep = "\t", header = FALSE, skip = 1, as.is = TRUE, stringsAsFactors = FALSE)
    col_num <- length(col_names)
    annoTable <- tmp[, 1:col_num]
    colnames(annoTable) <- col_names
    
    # filter out Gene_count != NA (Exon_count == NA)
    annoTable <- annoTable[!is.na(annoTable[, "Exon_count"]), , drop = FALSE]
    maf <- cbind(
        'Hugo_Symbol' = as.character(annoTable[,  "Gene_name"]),
        'Entrez_Gene_Id' = '1',
        'Center' = "GM",
        'NCBI_Build' = "hg19",
        'Chromosome' = as.character(annoTable[, "SV_chrom"]),
        'Start_position' = as.numeric(annoTable[, "SV_start"]),
        'End_position' = as.numeric(annoTable[, "SV_end"]),
        'Strand' = "+",
        'Variant_Classification' = apply(annoTable, 1, getSVC),
        'Variant_Type' = as.character(annoTable[, "SV_type"]),
        'Reference_Allele' = as.character(annoTable[, "REF"]),
        'Tumor_Seq_Allele1' = as.character(annoTable[, "ALT"]),
        'Tumor_Seq_Allele2' = "",                          
        'Tumor_Sample_Barcode' = sapply(annoTable[, "Samples_ID"], function(x) unlist(strsplit(x, ","))[2]),
        'Matched_Norm_Sample_Barcode' = sapply(annoTable[, "Samples_ID"], function(x) unlist(strsplit(x, ","))[1]),
        'Mutation_Status' = "Somatic",
        'VAF' = apply(annoTable, 1, getSVAF),
        'Frameshift' = as.character(annoTable[, "Frameshift"])
    )
    
    mafTemp <- as.data.frame(maf)
    return(mafTemp)
}

# get structural variant classification
getSVC <- function(mann){
    chrStart2 <- unlist(strsplit(mann["ALT"], "\\[|\\]"))[2]
    chr2 <- unlist(strsplit(chrStart2, ":"))[1]
    chr1 <- as.character(mann["SV_chrom"])
    if(mann["SV_type"] == "BND" && chr1 == chr2){
        return("Domestic_Ins")
    }else if(mann["SV_type"] == "BND" && chr1 != chr2){
        return("Translocation")
    }else if(mann["SV_type"] == "INV"){
        return("Inversion")
    }else if(mann["SV_type"] == "INS"){
        return("Foreign_Ins")
    }else if(mann["SV_type"] == "DUP"){
        return("Duplication")
    }else if(mann["SV_type"] == "DEL"){
        return("Deletion")
    }else{
        return("Unknown")
    }
}

# get structural variant allele frequency
getSVAF <- function(mann){
    normalSampleName <- as.character(unlist(strsplit(mann["Samples_ID"], ","))[1])
    tumorSampleName <- as.character(unlist(strsplit(mann["Samples_ID"], ","))[2])
    normalCounts <- unlist(strsplit(mann[normalSampleName], ":"))
    tumorCounts <- unlist(strsplit(mann[tumorSampleName], ":"))
    if (length(normalCounts) == 2) {
        n_counts <- as.numeric(unlist(strsplit(normalCounts[1], ",")))
        n_ref_count <- n_counts[1]
        n_alt_count <- n_counts[2]
        t_counts <- as.numeric(unlist(strsplit(tumorCounts[1], ",")))
        t_ref_count <- t_counts[1]
        t_alt_count <- t_counts[2]
        n_counts_2 <- as.numeric(unlist(strsplit(normalCounts[2], ",")))
        n_ref_count_2 <- n_counts_2[1]
        n_alt_count_2 <- n_counts_2[2]
        t_counts_2 <- as.numeric(unlist(strsplit(tumorCounts[2], ",")))
        t_ref_count_2 <- t_counts_2[1]
        t_alt_count_2 <- t_counts_2[2]
        n_vaf <- n_alt_count / (n_alt_count + n_ref_count)
        n_vaf_2 <- n_alt_count_2 / (n_alt_count_2 + n_ref_count_2)
        t_vaf <- t_alt_count / (t_alt_count + t_ref_count)
        t_vaf_2 <- t_alt_count_2 / (t_alt_count_2 + t_ref_count_2)
        VAF <- max(t_vaf, t_vaf_2)
    } else if(length(normalCounts) == 1) {
        n_counts <- as.numeric(unlist(strsplit(normalCounts[1], ",")))
        n_ref_count <- n_counts[1]
        n_alt_count <- n_counts[2]
        t_counts <- as.numeric(unlist(strsplit(tumorCounts[1], ",")))
        t_ref_count <- t_counts[1]
        t_alt_count <- t_counts[2]
        n_vaf <- n_alt_count / (n_alt_count + n_ref_count)
        t_vaf <- t_alt_count / (t_alt_count + t_ref_count)
        VAF <- t_vaf
    } else {
        VAF <- NULL
    }
    return(VAF)
}


#' Get CNV status from segmentation data
#'
#' @param df A dataframe - segmentation data
#' @param gene_list A string vector - list of human gene symbols.
#' @param return_num A logical - whether numeric value is returned or not. Default is \code{FALSE}.
#' @param return_text A logical value - if the return value shows text. Default value is \code{FALSE}.
#' @param reshape A logical value - if the result table is casted to matrix format. Default value is \code{FALSE}.
#' @param cutoff A numeric vector of length 2
#' @param geneInfoFile An imported geneInfo file (e.g., hg19, to replace hg18 in the original R package)
#'
#' @return A dataframe of CNV status
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang .data
#' @export

segment_to_cnv <- function(df,
                           gene_list,
                           return_num = FALSE,
                           return_text = FALSE,
                           reshape = FALSE,
                           cutoff = c(-0.3, 0.3),
                           geneInfoFile = NULL) {
    
    cnseg <- CNTools::CNSeg(df)
    
    if (is.null(geneInfoFile)) {
        utils::data("geneInfo", package = "CNTools")
    } else {
        geneInfo <- read.delim(geneInfoFile, sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
    }
    
    geneInfo <- geneInfo %>%
        filter(.data$genename %in% gene_list)
    
    rdByGene <- CNTools::getRS(cnseg,
                               by = "gene",
                               imput = FALSE,
                               XY = FALSE,
                               geneMap = geneInfo,
                               what = "median")
    
    cnv_res <- rdByGene@rs %>%
        `[`(6:ncol(rdByGene@rs)) %>%
        t()
    
    colnames(cnv_res) <- rdByGene@rs$genename
    
    cnv_res <- reshape2::melt(cnv_res)
    
    colnames(cnv_res) <- c("Sample", "Gene", "CNV")
    
    if (return_num == FALSE) {
        cnv_res %<>%
            mutate(CNV = test_cnv(.data$CNV,
                                  return_text = return_text,
                                  cutoff = cutoff))
    }
    
    if (reshape) {
        
        cnv_res <- reshape2::dcast(data = cnv_res,
                                   formula = Gene ~ Sample,
                                   value.var = "CNV")
        cnv_res <- tibble::column_to_rownames(cnv_res, var = "Gene")
    }
    
    if (is.null(geneInfoFile)) {
        rm(geneInfo, envir = .GlobalEnv)
    }
    return(cnv_res)
}

# Tools for processing CNV data

#' Determine if the segment mean indicates an amplification or a deletion
#'
#' @param values A numeric vector - the log2Ratio values for copy number changes.
#' @param cutoff A numeric vector of length 2 - the cutoff for calling amplifications or deletions
#' @param return_text A logical value - if the return value shows text. Default value is \code{FALSE}.
#'
#' @return A numeric vector of the same length with \code{values}
#' @export

test_cnv <- function(values,
                     cutoff = c(-0.3, 0.3),
                     return_text = FALSE) {
    
    if (cutoff[2] <= cutoff[1]) {
        stop("Cutoff values not accepted")
    }
    
    test_1 <- (values > cutoff[2])
    test_2 <- (values >= cutoff[1])
    
    test <- test_1 + test_2 - 1
    
    if (return_text) {
        
        test <- plyr::mapvalues(x = test,
                                from = c(1, 0, -1),
                                to = c("Gain", NA, "Loss"))
    }
    return(test)
}

filterMaf <- function(anned, tCount = 20, nCount = 10, tVaf = 0.02, nVaf = 0.02, minEvidence = 1, minTLOD, maxExactAll,
                      localRepeat = FALSE, dropOfftarget = FALSE, bed, genome = "BSgenome.Hsapiens.UCSC.hg19",
                      srCount = "(sample_ref_count)|(t_ref_count)", saCount = "(sample_alt_count)|(t_alt_count)",
                      rrCount = "(ref_ref_count)|(n_ref_count)", raCount = "(ref_alt_count)|(n_alt_count)",
                      sVafCol = "(tumor_f)|(VAF)") {
    
    srCountIndex <- grep(srCount, colnames(anned))
    saCountIndex <- grep(saCount, colnames(anned))
    if ((tCount != 0) && (length(srCountIndex) > 0) && (length(saCountIndex) > 0)) {
        anned <- anned[(as.numeric(anned[, srCountIndex[1]]) + as.numeric(anned[, saCountIndex[1]]))
                       > tCount, , drop = FALSE]
    }
    if(nrow(anned) == 0) return(anned)
    
    rrCountIndex <- grep(rrCount, colnames(anned))
    raCountIndex <- grep(raCount, colnames(anned))
    if ((nCount != 0) && (length(rrCountIndex) > 0) && (length(raCountIndex) > 0)) {
        anned <- anned[(as.numeric(anned[, rrCountIndex[1]]) + as.numeric(anned[, raCountIndex[1]]))
                       > nCount, , drop = FALSE]
    }
    if (nrow(anned) == 0) return(anned)
    
    if (tVaf != 0) {
        if (sVafCol %in% colnames(anned)){
            anned <- anned[as.numeric(anned[, sVafCol]) > tVaf, , drop = FALSE ]
        } else if ((length(srCountIndex) > 0) && (length(saCountIndex) > 0)) {
            t_f <- as.numeric(anned[, saCountIndex[1]]) / 
                (as.numeric(anned[, srCountIndex[1]]) + as.numeric(anned[, saCountIndex[1]]))
            anned <- anned[t_f > tVaf, , drop = FALSE]
        } else {
            
        }
    }
    if (nrow(anned) == 0) return(anned)
    
    if ((nVaf != 1) && (length(rrCountIndex) > 0) && (length(raCountIndex) > 0)) {
        normal_f <- as.numeric(anned[, raCountIndex[1]]) / 
            (as.numeric(anned[, rrCountIndex[1]]) + as.numeric(anned[, raCountIndex[1]]))
        anned <- anned[normal_f < nVaf, , drop = FALSE]
    }
    if(nrow(anned) == 0) return(anned)
    
    if((minEvidence > 1) & ("evidence" %in% colnames(anned))){ # Sean added
        evidenceNum <- lapply(anned[, "evidence"], function(x) length(unlist(strsplit(x, " \\| "))))
        evidenceNum <- unlist(evidenceNum)
        anned <- anned[evidenceNum >= minEvidence, , drop = FALSE]
    }
    
    if(!missing(minTLOD) & length(grep("t_lod_fstar", colnames(anned))) != 0){
        anned <- anned[as.numeric(anned[, "t_lod_fstar"]) >= minTLOD | is.na(anned[, "t_lod_fstar"] | anned[, "t_lod_fstar"] == 'NA' | anned[, "t_lod_fstar"] == '.' | anned[, "t_lod_fstar"] == ''), , drop = FALSE]
    }
    
    if(!missing(maxExactAll)){
        anned[anned[, "exac_all"] %in% c(NA, "NA", "", "."), "exac_all"] <- 0
        anned <- anned[as.numeric(anned[, "exac_all"]) < maxExactAll,  , drop = FALSE]
    }
    
    if(localRepeat){
        anned <- cbind(anned, getUpSeqs(anned, chrom = "Chromosome", end = "End_position", flanking = 25, genome = genome),
                       getDownSeqs(anned, chrom = "Chromosome", end = "End_position", flanking = 25, genome = genome))
        anned <- anned[anned[, "local_repeat_u"] != "TRUE" & anned[, "local_repeat_d"] != "TRUE", , drop = FALSE]
        # anned <- anned[anned[, "local_repeat_u"] != "TRUE", , drop = FALSE]
        # anned <- anned[anned[, "local_repeat_d"] != "TRUE", , drop = FALSE]
    }
    
    if(dropOfftarget && !missing(bed)){
        anned <- filterByBed(anned, bed, chrom = "Chromosome", start = "Start_position", end = "End_position")
    }
    
    return(anned)
}

filterByBed <- function(anned, bed, chrom = "chrom", start = "start", end = "end"){
    annedGR <- mat2GR(anned, chrom = chrom, start = start, end = end, what = "snp")
    bedGR <- mat2GR(bed, chrom = 1, start = 2, end = 3, what = "snp")
    ol <- findOverlaps(annedGR, bedGR, type = "any")
    anned <- anned[queryHits(ol), ]
    
    return(anned)	
}

mat2GR <- function(mat, chrom = "chrom", start = "start", 
                   end, strand, what = "snp"){
    require("GenomicRanges")
    exclu <- c(chrom, start)
    if(missing(strand) || is.null(strand) || !strand %in% colnames(mat)){
        stra <- factor(rep("+", nrow(mat)))
    }else {
        stra <- mat[[strand]]
        exlcu <- c(exclu, strand)
    }
    if(missing(end) || is.null(end)){
        end <- start
    }else{
        exclu <- c(exclu, end)
    }
    if(class(mat)[1] == "data.table"){
        if(length(grep("chr", mat[[chrom]])) == 0){
            chroms <- paste("chr", gsub(" ", "", as.character(mat[[chrom]])), sep = "")
        }
        argList <- list(seqnames = chroms, 
                        ranges = IRanges(start = mat[[start]], end = mat[[end]]),
                        strand = stra)
        for(mn in setdiff(colnames(mat), exclu)){
            argList[[tolower(mn)]] <- mat[[mn]]
        }
    }else{
        if(length(grep("chr", mat[, chrom])) == 0){
            mat[, chrom] <- paste("chr", gsub(" ", "", as.character(mat[, chrom])), sep = "")
        }
        argList <- list(seqnames = mat[, chrom], 
                        ranges = IRanges(start = as.numeric(mat[, start]), 
                                         end = as.numeric(mat[, end])), strand = stra)
        for(mn in setdiff(colnames(mat), exclu)){
            argList[[tolower(mn)]] <- mat[, mn]
        }
    }
    argList[["dataType"]] <- rep(what, nrow(mat))
    gr <- do.call(GRanges, args = argList)
    
    return(gr)
}

### gets the somatic or germline pindels
getPindelSP <- function(anned, what = c("somatic", "germline"), sName, isFile = TRUE){
    what = match.arg(what)
    if(isFile){
        annedMat <- read.delim(anned, sep = "\t", header = TRUE, as.is = TRUE)
    }else{
        annedMat <- anned
    }
    colnames(annedMat) <- tolower(colnames(annedMat))
    if(missing(sName)){
        sName <- paste(unique(annedMat[, "sample_name"]), sep = "", collapse = "|")
        if (length(sName) == 0) sName <- NA
    }
    if(what == "germline"){
        return(annedMat[-grep("NA;|;NA", annedMat[, "samples"]),])
    }
    annedMat <- annedMat[as.numeric(annedMat[, "nsample"]) == 1, , drop=F]
    pos <- grep(sName, annedMat[, "samples"])
    if(is.na(pos) || length(pos) == 0){
        print("No somatic variants found")
        return(NULL)
    }
    return(annedMat[pos, , drop=F])
}

getUpSeqs <- function (mat, chrom = "chrom", end = "end", flanking = 25,
                       genome = "BSgenome.Hsapiens.UCSC.hg19", threshold = 4) {
    flanked <- getFlankingSeqs(mat, chrom = chrom, start = end, end = end, what = "S",
                               flanking = flanking, genome = genome)
    upseqs <- cbind(up_seqs = flanked[, "left"],
                    local_repeat_u = sapply(sapply(flanked[, "left"],
                                                   function(x) paste(rev(unlist(strsplit(x, "*"))), collapse = "")),
                                            hasRepeats, threshold = threshold))
    return(upseqs)
}

getDownSeqs <- function (mat, chrom = "chrom", end = "end", flanking = 25,
                         genome = "BSgenome.Hsapiens.UCSC.hg19", threshold = 4) {
    flanked <- getFlankingSeqs(mat, chrom = chrom, start = end, end = end, what = "S",
                               flanking = flanking, genome = genome)
    downseqs <- cbind(down_seqs = flanked[, "right"],
                      local_repeat_d = sapply(flanked[, "right"],
                                              hasRepeats, threshold = threshold))
    return(downseqs)
}

getFlankingSeqs <- function(mutMat, flanking = 20, what = c("S", "D", "I"), 
                            chrom = "chrom", start = "start", end = "end", 
                            genome = "BSgenome.Hsapiens.UCSC.hg19"){
    require(Biostrings)
    what <- match.arg(what)
    require(genome, character.only = TRUE)
    seqs <- get(genome, pos = grep(genome, search()))
    if(length(grep("chr", mutMat[, chrom])) == 0){
        mutMat[, chrom] <- paste("chr", mutMat[, chrom], sep = "")
    }
    mutMat[mutMat[, chrom] == "chrMT", chrom] <- "chrM"
    if(flanking != 0){
        if(what %in% c("D", "S")){
            left <- mat2GR(cbind(chrom = gsub(" ", "", mutMat[, chrom]), 
                                 start = as.numeric(mutMat[, start]) - flanking, 
                                 end = as.numeric(mutMat[, start]) - 1), chrom = "chrom", start = "start", end = "end")
            right <- mat2GR(cbind(chrom = gsub(" ", "", mutMat[, chrom]), start = as.numeric(mutMat[, end]) + 1, 
                                  end = as.numeric(mutMat[, end]) + flanking), chrom = "chrom", start = "start", end = "end")
            
        }else{
            left <- mat2GR(cbind(chrom = gsub(" ", "", mutMat[, chrom]), 
                                 start = as.numeric(mutMat[, start]) - flanking + 1, 
                                 end = as.numeric(mutMat[, start]) ), chrom = "chrom", start = "start", end = "end")
            right <- mat2GR(cbind(chrom = gsub(" ", "", mutMat[, chrom]), 
                                  start = as.numeric(mutMat[, end]), 
                                  end = as.numeric(mutMat[, end]) + flanking - 1), chrom = "chrom", start = "start", end = "end")
        }
        return(cbind(left = Biostrings:::getSeq(seqs, left, as.character = TRUE), 
                     right = Biostrings:::getSeq(seqs, right, as.character = TRUE)))
    }
    return(Biostrings:::getSeq(seqs, mat2GR(cbind(chrom = gsub(" ", "", mutMat[, chrom]), 
                                                  start = as.numeric(mutMat[, start]), end = as.numeric(mutMat[, end])), chrom = "chrom", 
                                            start = "start", end = "end"), as.character = TRUE))									 	
}

hasRepeats <- function(str, threshold = 4){
    pos <- gregexpr(substr(str, 1, 1), str)
    
    if(length(pos[[1]]) == 1){ return(FALSE) }
    
    rep <- substr(str, as.numeric(pos[[1]][1]), as.numeric(pos[[1]][2]) - 1)
    repPos <- gregexpr(rep, str, perl = TRUE)
    
    if(length(repPos[[1]]) >= threshold){
        pWidth <- as.numeric(repPos[[1]])[2:length(repPos[[1]])] - as.numeric(repPos[[1]])[1:(length(repPos[[1]]) - 1)]
        if(all(1:(threshold-1) %in% grep(paste("^", pWidth[1], "$", sep = ""), pWidth))){
            return(TRUE)
        }
    }
    return(FALSE)
}

# Construct copy-number (CN) spectra plot
cnSpecPlot_byChrome <- function(cnSeg, isFile = TRUE, chrome = 'all', cutoff = FALSE, organism = 'human', genome = 'hg19') {
    if (is.null(chrome)) {
        stop("No chrome selected!")
    }
    require("GenVisR")
    
    # chrome information
    if (organism == "human") {
        require("BSgenome.Hsapiens.UCSC.hg19")
        genome.lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
        genome = "hg19"
    } else if (organism == "mouse") {
        require("BSgenome.Mmusculus.UCSC.mm10")
        genome.lengths <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
        genome = "mm10"
    } else {
        stop("organism not supported yet!!!")
    }
    
    # reading CN segments
    if (isFile) {
        cnSeg_input <- read.delim(cnSeg, sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
        cnSeg_GenVisR <- as.data.frame(cnSeg_input[, c("chrom", "loc.start", "loc.end", "seg.mean", "ID")])
        colnames(cnSeg_GenVisR) <- c("chromosome", "start", "end", "segmean", "sample")
    } else {
        cnSeg_GenVisR <- cnSeg
    }
    
    # CN segment filtering
    if (length(cutoff) == 1) {
        if (cutoff == TRUE) {
            # Default settings
            cnSeg_GenVisR <-  cnSeg_GenVisR[(as.numeric(cnSeg_GenVisR[, "segmean"]) > 0.3 |
                                                 as.numeric(cnSeg_GenVisR[, "segmean"]) < -0.3), ]
        }
    } else if (length(cutoff) == 2) {
        if (as.numeric(cutoff[1]) <= as.numeric(cutoff[2])) {
            cnSeg_GenVisR <-  cnSeg_GenVisR[(as.numeric(cnSeg_GenVisR[, "segmean"]) > cutoff[2] |
                                                 as.numeric(cnSeg_GenVisR[, "segmean"]) < cutoff[1]), ]
        } else {
            stop("Cutoff values not accepted!")
        }
    }
    
    # CN log ratio scaling
    seg.mean <- as.numeric(cnSeg_GenVisR[, "segmean"])
    cnSeg_GenVisR[, "segmean"] <- if_else(abs(seg.mean) >= 2, 2 * sign(seg.mean), seg.mean)
    
    # color legends (can be changed as settings later!!!)
    plotLayer <- scale_fill_gradientn(colors = c("blue", "white", "red"),
                                      values = scales::rescale(c(-2, 0, 2)),
                                      limits = c(-2, 2),
                                      name = "log_ratio")
    
    # different chrome combination allowed
    if (length(chrome) == 1 && chrome == 'all') {
        cnSpec(cnSeg_GenVisR, genome = genome, plotLayer = plotLayer, CNscale = "relative")
    } else {
        cnSeg_GenVisR.sub <- NULL
        for (chr in chrome) {
            cnSeg_GenVisR.sub <- rbind(cnSeg_GenVisR.sub, cnSeg_GenVisR[cnSeg_GenVisR[, "chromosome"] == chr, ])
        }
        genome.range.sub <- as.data.frame(list("chromosome" = chrome,
                                               "start" = rep('0', length(chrome)),
                                               "end" = genome.lengths[paste0('chr', chrome)]))
        cnSpec(cnSeg_GenVisR.sub, y = genome.range.sub, plotLayer = plotLayer, genome = genome, CNscale = "relative")
    }
}

# Construct copy-number (CN) frequency plot
cnFreqPlot_byChrome <- function(cnvFile, chrome = 'all', cutoff = FALSE, organism = 'human', genome = 'hg19') {
    if (is.null(chrome)) {
        stop("No chrome selected!")
    }
    require("GenVisR")
    
    # chrome information
    if (organism == "human") {
        genome = "hg19"
    } else if (organism == "mouse") {
        genome = "mm10"
    } else {
        stop("organism not supported yet!!!")
    }
    
    # reading CN segments
    cnSeg_input <- read.delim(cnvFile, sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
    cnSeg_GenVisR <- as.data.frame(cnSeg_input[, c("chrom", "loc.start", "loc.end", "seg.mean", "ID")])
    colnames(cnSeg_GenVisR) <- c("chromosome", "start", "end", "segmean", "sample")
    
    # CN segment filtering
    if (length(cutoff) == 1) {
        if (cutoff == TRUE) {
            # Default settings
            cnSeg_GenVisR <-  cnSeg_GenVisR[(as.numeric(cnSeg_GenVisR[, "segmean"]) > 0.3 |
                                                 as.numeric(cnSeg_GenVisR[, "segmean"]) < -0.3), ]
        }
    } else if (length(cutoff) == 2) {
        if (as.numeric(cutoff[1]) <= as.numeric(cutoff[2])) {
            cnSeg_GenVisR <-  cnSeg_GenVisR[(as.numeric(cnSeg_GenVisR[, "segmean"]) > cutoff[2] |
                                                 as.numeric(cnSeg_GenVisR[, "segmean"]) < cutoff[1]), ]
        } else {
            stop("Cutoff values not accepted!")
        }
    }
    
    # CN log ratio scaling
    seg.mean <- as.numeric(cnSeg_GenVisR[, "segmean"])
    cnSeg_GenVisR[, "segmean"] <- if_else(abs(seg.mean) >= 2, 2 * sign(seg.mean), seg.mean)
    
    # convert to absolute values
    cnSeg_GenVisR[, "segmean"] <- cnSeg_GenVisR[, "segmean"] + 2
    
    # different chrome combination allowed
    if (length(chrome) == 1 && chrome == 'all') {
        cnFreq(cnSeg_GenVisR, plotType = "frequency", genome = genome)
    } else {
        cnSeg_GenVisR.sub <- NULL
        for (chr in chrome) {
            cnSeg_GenVisR.sub <- rbind(cnSeg_GenVisR.sub, cnSeg_GenVisR[cnSeg_GenVisR[, "chromosome"] == chr, ])
        }
        cnFreq(cnSeg_GenVisR.sub, plotChr = paste0('chr', chrome), plotType = "frequency", genome = genome)
    }
}

# Format (parse, dedup and sort) chrome list from a string, e.g., "1-3, 4-7, 8-12, 13-19, 20-22,X,Y"
chromeList <- function(chrome, chromeListFull = c(1:22, 'X', 'Y'), organism = 'human') {
    # chrome information
    if (organism == "human") {
        chromeListFull = c(1:22, 'X', 'Y')
    } else if (organism == "mouse") {
        chromeListFull = c(1:19, 'X', 'Y')
    } else {
        stop("organism not supported yet!!!")
    }
    chrome <- gsub(' ', '', chrome)
    chrome <- gsub('chr', '', chrome, ignore.case = TRUE)
    if (sum(grepl('all', chrome, ignore.case = TRUE)) > 0) {
        return('all')
    } else {
        strs <- unlist(strsplit(chrome, ',|;'))
        chrList <- NULL
        for (chr in strs) {
            if (grepl('-|:|~', chr)) {
                chrs <- unlist(strsplit(chr, '-|:|~'))
                chr <- chrs[1]:chrs[2]
            }
            chrList <- union(chrList, chr)
        }
        chrList <- intersect(chromeListFull, chrList)
        if (length(chrList) == length(chromeListFull)) {
            return('all')
        } else if (length(chrList) == 0) {
            return(NULL)
        } else {
            return(chrList)
        }
    }
}

readColorCode <- function(feature_color_file) {
    color_table <- read.delim(feature_color_file, sep = ",", header = TRUE, as.is = TRUE)
    
    feature_list <- unique(color_table[, "Feature"])
    feature_color <- list()
    for (feature in feature_list) {
        class_list <- color_table[color_table[, "Feature"] == feature, "Class"]
        color_list <- color_table[color_table[, "Feature"] == feature, "Color"]
        names(color_list) = class_list
        feature_color <- c(feature_color, list(color_list))
    }
    names(feature_color) <- feature_list
    
    return(feature_color)
}
