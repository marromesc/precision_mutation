# Maria Roman Escorza - 2023 01 11

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

library(reshape2)
library(ggplot2)

meta_datapath <- './results/SampleSheet.csv'
SampleSheet <- read.csv(meta_datapath)

summaryFigure <- function(x, col=NULL){
  x$status <- ifelse(x$case_control =='ipsilateral DCIS', 'DCIS-Case', ifelse(x$case_control=='case', 'INV-Case', 'Control'))

  x$sample_origin <- paste0(x$cohort, '-', x$platform)
  x$her2 <- factor(ifelse(x$her2 == 0, "Negative", ifelse(x$her2 == 1, "Positive", "Unknown")), levels = c("Positive","Negative", 'Unknown'), labels=c("HER2+","HER2-", 'U'))
  x$er <- factor(ifelse(x$er == 0, "Negative", ifelse(x$er == 1, "Positive", "Unknown")), levels=c("Positive","Negative", 'Unknown'),labels=c("ER+","ER-", 'U'))
  x$grade <- factor(ifelse(x$grade %in% c("999", 'NoData') | is.na(x$grade), 'Unknown', x$grade),levels=c('3','2', '1','Unknown'), labels=c('3','2', '1','U'))
  x$age_diagnosis <- factor(ifelse(x$age_diagnosis==999, "Unknown", ifelse(x$age_diagnosis>50, "Yes", ifelse(x$age_diagnosis<50, "No", "Unknown"))), levels=c("Yes","No","Unknown"), labels=c(">50", "<50", 'U'))
  x$radiotherapy <- factor(ifelse(x$radiotherapy==999, "Unknwon", ifelse(x$radiotherapy==1, "RT", ifelse(x$radiotherapy==0, "No-RT", "Unknown"))),levels=c('RT', "No-RT", 'Unknown'), labels=c('RT', "No-RT", 'U'))
  
  if(!is.null(col)){
    for(i in col){
      x[,i] <- ifelse(is.na(x[,i]) | x[,i] == 'ICNA', 'Unknown', x[,i])
      x[,i] <- factor(x[,i], 
                      levels = unique(ifelse(is.na(x[,i]), 'Unknown', x[,i])), 
                      labels = ifelse(unique(x[,i]) == 'LumA', 'A', ifelse(unique(x[,i]) == 'LumB', 'B', ifelse(unique(x[,i]) == 'Her2', 'H', ifelse(unique(x[,i]) == 'Basal', 'Ba', ifelse(unique(x[,i]) == 'Her2', 'H', ifelse(unique(x[,i]) == 'Normal', 'N', ifelse( unique(x[,i]) == 'Unknown' , 'U', unique(x[,i]))))))))
      )
    }
    e <-x[,c('sample_origin', "her2","er","grade","age_diagnosis","radiotherapy",col, "status")]
  } else {
    e<-x[,c('sample_origin',"her2","er","grade","age_diagnosis","radiotherapy", "status")]
  }
  
  colnames(e)[ncol(e)] <- 'Class'
  
  df<-character()
  for (g in colnames(e)){
    message(g)
    df <- rbind(df,cbind(melt(prop.table(table(e$Class,e[,g]),margin = 1)),g))
  }
  df<-df[complete.cases(df),]
  colnames(df)<-c("Class","value","prop","variable")
  
  df<-df[df$prop!=0,]
  df<-df[df$variable!="Class",]
  
  df$value<-factor(df$value,levels = rev(unique(df$value)))
  
  df$type <- "ALL"
  #df[df$variable %in% c('cn', 'rna', 'mutation'),'type'] <- 'Data type'
  df[df$variable %in% c("sample_origin"),"type"]<-"Cohort"
  df[df$variable %in% c("her2",'er', "grade"),"type"]<-"Histological"
  df[df$variable %in% c("age_diagnose"),"type"]<-"Clinical"
  df[df$variable %in% c("radiotherapy"),"type"]<-"Treatment"
  #df[df$variable %in% col,"type"]<-"Subtype"
  df$type<-factor(df$type,levels=c('Cohort',"Histological", 'Clinical', "Treatment"))
  df<-df[rev(order(df$value)),]
  
  cols.data.type <- c('#393ecf', '#a0a2e8')
  cols.cohort <- c('#d53e4f','#fc8d59', '#fee08b', '#e6f598', '#99d594', '#3288bd')
  cols.er        <- c("#74c476", "#bae4b3")  
  cols.her2      <- c("#74c476", "#bae4b3")
  cols.grade     <- c("#33cec8","#91e4e2", "#f0fbfb")
  cols.age <- c('lightpink','lightcoral')
  cols.treatment <- c('#fcbba1','#fc9272')
  
  vect <- c(cols.cohort,cols.er,cols.her2,cols.grade,cols.treatment,cols.age)
  if(!is.null(col)){
    for (i in col){
      vect <- c(vect, RColorBrewer::brewer.pal(length(unique(x[,i])), 'Set3'))
    }
  } 
  cols <- rev(vect)
  
  df$variable<-factor(df$variable,levels=rev(c("sample_origin", "er","her2", "grade","age_diagnose",'radiotherapy',col)))
  
  cols <- c('white', cols)
  names(cols)<-c('U', levels(df$value)[-which(levels(df$value) == 'U')])
  
  df$label <- as.character(df$value)
  df$label <- gsub("^ER\\+|^HER2\\+|^PR\\+","\\+",df$label)
  df$label <- gsub("^ER\\-|^HER2\\-|^PR\\-","\\-",df$label)
  #df$label <- ifelse(df$label == 'KCL-SNParray', 'K1', ifelse(df$label == 'KCL-lpWGS', 'K2', ifelse(df$label == 'NKI-lpWGS', 'N1', df$label)))
  
  df$Class <- ifelse(df$Class == 'Control', 
                     paste0('Control (n=',as.numeric(table(x$status)[names(table(x$status)) == 'Control']),')'), 
                     ifelse(df$Class == 'INV-Case', 
                            paste0('INV-Case (n=', as.numeric(table(x$status)[names(table(x$status)) == 'INV-Case']),')'), 
                            paste0('DCIS-Case (n=', as.numeric(table(x$status)[names(table(x$status)) == 'DCIS-Case']),')')))
  
  sum <- ggplot(df,aes(y=prop,x=variable))+
    geom_bar(stat="identity", position="fill",aes(fill=value),color="gray30",size=0.2)+
    scale_fill_manual(values = cols)+
    geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
              position = position_stack(vjust = 0.5)) +
    facet_grid(type~Class, scales="free", space = "free")+
    labs(y="Proportion of cases")+
    scale_x_discrete(expand = c(0, 0), breaks= levels(df$variable),
                     labels= rev(
                       c( "Cohort", "ER","Her2", "Grade", 'Age of diagnosis', "Radiotherapy",col)))+
    scale_y_continuous(breaks=c(0,0.5,1),labels = c("0%","50%","100%"))+
    theme_grey(base_size = 13) + 
    theme(
      axis.text=element_text(size = 13*0.7, face="bold",colour = "black"),
      axis.title=element_blank(),
      axis.ticks=element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      panel.background = element_blank(),
      strip.text.x = element_text(size = 13*0.8, face="bold"),
      strip.text.y = element_blank(),
      panel.spacing.x = unit(0.00001, "lines"),
      panel.spacing.y = unit(0.8, "lines"))+
    coord_flip()
  
  return(sum)
}

#main figure
pdf('./results/summary_samples.pdf', height = 5, width = 11)
print(summaryFigure(SampleSheet))
dev.off()

