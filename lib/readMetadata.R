readMetadata <- function(datapath){
  metadata1 <- as.data.frame(readxl::read_xlsx(datapath))
  metadata1 <- metadata1[metadata1$first_subseq_event %in% c('ipsilateral DCIS', 'ipsilateral IBC', 'death', 'NA'),]
  metadata2 <- as.data.frame(readxl::read_xlsx(datapath, sheet = 2))
  metadata2 <- metadata2[metadata2$redcap_event == 'PRI',]
  metadata_merged <- as.data.frame(merge(metadata1, metadata2))
  
  return(metadata_merged)
}