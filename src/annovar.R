library(foreach)
require(doMC)
registerDoMC(40)

setwd('/mnt/albyn/maria/precision_mutation')

# KCL Panel ---------------------------------------------------------------

system('mkdir ./data/Panel/DCIS_Precision_Panel_KCL')

path_in <- "/home/maria/albyn/prj_precision/targetedseq/Sloane_Vandna/VCF/PDCIS/"
vcfFiles <- list.files(path=path_in, pattern=paste("*_P_bcf_processed.vcf"), full.names = FALSE)

foreach (file=vcfFiles) %dopar%
  {
    message(file)
    Command <- "perl /mnt/albyn/common/annovar/table_annovar.pl --buildver hg19 -vcfinput "
    Command <- paste0(Command, path_in, file, ' -protocol refGene,exac03,cosmic70,avsnp147,dbnsfp30a,esp6500siv2_all,ALL.sites.2015_08,clinvar_20170905 -operation g,f,f,f,f,f,f,f -polish /mnt/albyn/common/annovar/humandb')
    system(Command)
  }

Command <- 'ln -s '
Command <- paste0(Command, path_in, '*_bcf_processed.vcf.hg19_multianno.txt ', getwd(), '/data/Panel/DCIS_Precision_Panel_KCL')
system(Command)


# NKI Panel ---------------------------------------------------------------

system('mkdir ./data/Panel/DCIS_Precision_Panel_NKI')

path_in <- "/home/maria/albyn/prj_precision/targetedseq/PRECISION_NKI/ANNO_VCF/"
vcfFiles <- list.files(path=path_in, pattern=paste("*.vcf"), full.names = FALSE)

foreach (file=vcfFiles) %dopar%
  {
    message(file)
    id <- gsub('.vcf', '', file)
    
    Command <- "perl /mnt/albyn/common/annovar/table_annovar.pl --buildver hg19 -vcfinput "
    Command <- paste0(Command, path_in, file, ' -protocol refGene,exac03,cosmic70,avsnp147,dbnsfp30a,esp6500siv2_all,ALL.sites.2015_08,clinvar_20170905 -operation g,f,f,f,f,f,f,f -polish /mnt/albyn/common/annovar/humandb')
    system(Command)
  }

Command <- 'ln -s '
Command <- paste0(Command, path_in, '*.vcf.hg19_multianno.txt ', getwd(), '/data/Panel/DCIS_Precision_Panel_NKI')
system(Command)
