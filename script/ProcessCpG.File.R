setwd("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/New_analysis_Bowtie2/Sorted_Bam/")
my.filelist<-list.files(pattern = ".CpG.txt")
my.header<-unlist(lapply(strsplit(my.filelist,"_"),'[[',1))
for (i in 1:length(my.filelist)){
  x<-read.table(my.filelist[i],header = T)
  x$chrBase<-paste0("chr",x$chrBase)
  x$chr<-paste0("chr",x$chr)
  xx<-x[grep("^chr[1-9,X,Y,MT]",x$chr),]
  write.table(xx,
              file = paste0(my.header[i],"_CpGs.chr.txt"),
              quote = F,
              row.names = F,
              sep = "\t")
  print(i)
}

