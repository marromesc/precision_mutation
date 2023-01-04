readMutect <- function(datapath, source, pattern = NULL){
  # combine all mutects 
  if (is.null(pattern)){
    ls.files = list.files(path = c(datapath), pattern = '.tsv', full.names = TRUE)
  } else {
    ls.files = do.call(rbind, lapply(pattern, function(x) list.files(path = c(datapath), pattern = x, full.names = TRUE)))[,1]
  }
  
  datlist = list()
  for(i in 1:length(ls.files)){
    #message(ls.files[i])
    if( file.info(ls.files[i])$size > 0 ){
      dat.mat <- read.delim2(file = ls.files[i]) %>% 
        tbl_df() %>% 
        mutate(tool = "mutect", source = source)
      
      datlist[[i]] = dat.mat
    } 
  }
  
  df_comb = do.call(rbind, datlist);dim(df_comb) 
  return(df_comb)
}