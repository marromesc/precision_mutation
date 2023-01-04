readMetadata <- function(datapath){
  metadata1 <- as.data.frame(readxl::read_xlsx(datapath))
  metadata2 <- as.data.frame(readxl::read_xlsx(datapath, sheet = 2))
  metadata2 <- metadata2[metadata2$redcap_event == 'PRI',]
  metadata_merged <- as.data.frame(merge(metadata1, metadata2))
  
  stopifnot(nrow(metadata2) == nrow(metadata_merged))
  return(metadata_merged)
}