# Maria Roman Escorza - 2023 02 08

# Load libraries and data path --------------------------------------------

setwd('/mnt/albyn/maria/precision_mutation')

system('mkdir ./results/clinical_association')

mutation_datapath <- './results/Filtered_Mutations_Compiled_5%.csv'
meta_datapath <- './results/SampleSheet.csv'

eventDataFrame <- read.csv(mutation_datapath)
eventDataFrame <- eventDataFrame[!is.na(eventDataFrame$case_control),]
SampleSheet <- read.csv(meta_datapath)

GenesPanel <- readRDS('./data/GenesPanel.RDS')


# Gene Matrix -------------------------------------------------------------

# extract patients
patients <- SampleSheet$patient_id
  
# prepare matrix
geneMatrix <- matrix("",nrow=length(patients),ncol=length(GenesPanel))
rownames(geneMatrix) <- patients
colnames(geneMatrix) <- GenesPanel
ER <- rep("",length(patients));  Her2 <- rep("",length(patients));  Grade <- rep("",length(patients)); RT = rep('', length(patients)); Batch = rep('', length(patients))

for (i in 1:nrow(geneMatrix)) {
  Ind <- which(eventDataFrame$Entrez_Gene_Id==rownames(geneMatrix)[i])
  pat_id <- rownames(geneMatrix)[i]
  ER[i] <- SampleSheet$er[which(SampleSheet$patient_id == pat_id)]
  Her2[i] <- SampleSheet$her2[which(SampleSheet$patient_id == pat_id)]
  Grade[i] <- SampleSheet$grade[which(SampleSheet$patient_id == pat_id)]
  RT[i] <- SampleSheet$radiotherapy[which(SampleSheet$patient_id == pat_id)]
  Batch[i] <-  SampleSheet$cohort[which(SampleSheet$patient_id == pat_id)]
  for (j in 1:ncol(geneMatrix)) {
    Ind <- which(eventDataFrame$Entrez_Gene_Id==rownames(geneMatrix)[i] & eventDataFrame$Hugo_Symbol==colnames(geneMatrix)[j])
    if (length(Ind)==1) {
      geneMatrix[i,j] <- eventDataFrame$Consequence[Ind]
    }
    if (length(Ind)>1) {
      if (length(unique(eventDataFrame$Entrez_Gene_Id[Ind]))==1) {  
        geneMatrix[i,j] <- "multi_hit"
      }
      else {
          geneMatrix[i,j] <- eventDataFrame$Consequence[Ind[1]]
      }
  }
  }
}

Grade <- ifelse(Grade == 3, 'High', ifelse(Grade == 2, 'Int', ifelse(Grade == 1, 'Low', 'Unknown')))
ER <- ifelse(ER == 1, 'Positive', ifelse(ER == 0, 'Negative', 'Unknown'))
RT <- ifelse(RT == 1, 'RT+', ifelse(RT == 0, 'RT-', 'Unknown'))
Her2 <- ifelse(Her2 == 1, 'Positive', ifelse(Her2 == 0, 'Negative', 'Unknown'))


# Clinical association ----------------------------------------------------

write.table(data.frame(Gene="gene",what="what", Positive="Positive",Positive_tot="Positive_tot", Negative="Negative",Negative_Tot="Negative_Tot",Pval="p.value"), 
            file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=FALSE)
for (g in GenesPanel) {
  Res <- data.frame(Gene=g,what="ER", Positive=sum(geneMatrix[which(ER=="Positive"),g]!=""),Positive_tot=length(which(ER=="Positive")), Negative=sum(geneMatrix[which(ER=="Negative"),g]!=""),Negative_Tot=length(which(ER=="Negative")))
  Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
  write.table(Res, file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
  Res <- data.frame(Gene=g,what="Her2", Positive=sum(geneMatrix[which(Her2=="Positive"),g]!=""),Positive_tot=length(which(Her2=="Positive")), Negative=sum(geneMatrix[which(Her2=="Negative"),g]!=""),Negative_Tot=length(which(Her2=="Negative")))
  Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
  write.table(Res, file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
}
Res <- data.frame(Gene="all genes",what="ER", Positive=sum(geneMatrix[which(ER=="Positive"),]!=""),Positive_tot=length(GenesPanel)*length(which(ER=="Positive")), Negative=sum(geneMatrix[which(ER=="Negative"),]!=""),Negative_Tot=length(GenesPanel)*length(which(ER=="Negative")))
Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
write.table(Res, file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
Res <- data.frame(Gene="all genes",what="Her2", Positive=sum(geneMatrix[which(Her2=="Positive"),]!=""),Positive_tot=length(GenesPanel)*length(which(Her2=="Positive")), Negative=sum(geneMatrix[which(Her2=="Negative"),]!=""),Negative_Tot=length(GenesPanel)*length(which(Her2=="Negative")))
Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
write.table(Res, file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)

write.table(data.frame(Gene="gene",what="what", Positive="Int/Low",Positive_tot="Int/Low_tot", Negative="High",Negative_Tot="High_Tot",Pval="p.value"), 
            file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
for (g in GenesPanel) {
  Res <- data.frame(Gene=g,what="Grade", Positive=sum(geneMatrix[which(Grade %in% c("Low","Intermediate","1","2")),g]!=""),Positive_tot=length(which(Grade %in% c("Low","Intermediate","1","2"))), Negative=sum(geneMatrix[which(Grade %in% c("High","3")),g]!=""),Negative_Tot=length(which(Grade %in% c("High","3"))))
  Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
  write.table(Res, file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
}
Res <- data.frame(Gene="all genes",what="Grade", Positive=sum(geneMatrix[which(Grade %in% c("Low","Intermediate","1","2")),]!=""),Positive_tot=length(GenesPanel)*length(which(Grade%in% c("Low","Intermediate","1","2"))), Negative=sum(geneMatrix[which(Grade %in% c("High","3")),]!=""),Negative_Tot=length(GenesPanel)*length(which(Grade %in% c("High","3"))))
Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
write.table(Res, file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)

# correct fdr
Res <- read_tsv("./results/clinical_association/ClinicalMutationsStatistics.csv")
Res$fdr <- p.adjust(as.numeric(Res[['p.value']]), method = 'fdr')
write.table(Res, file="./results/clinical_association/ClinicalMutationsStatistics.csv", row.names=FALSE,sep="\t", quote = FALSE)
