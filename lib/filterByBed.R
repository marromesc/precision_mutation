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