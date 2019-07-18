if (!require(DSS)){
  BiocManager::install("DSS")
}
require(bsseq)
options(stringsAsFactors = F)
setwd("/home/cyz/Bigstore/BigData/analysis_work/Amer/Amer_run/Sorted_BAM/coverage_file/")
# prepressing file, they contain alternative chromosome
# my.filelist<-list.files(pattern = "cov.gz")
# my.filename<-unlist(lapply(strsplit(my.filelist,"_"),"[",1))
# for (i in (1:length(my.filelist))){
#   #read in file
#   tmp<-read.table(my.filelist[i])
#   # filter out alternative chromosome
#   tmp.1<-tmp[grep("^[1-9,X,Y,M]",tmp$V1),]
#   # add chrX instead of X
#   tmp.1[,1]<-paste0("chr",tmp.1$V1)
#   tmp.2<-cbind.data.frame(chr=tmp.1$V1,pos=tmp.1$V2,N=c(tmp.1$V5+tmp.1$V6),X=tmp.1$V5)
#   write.table(tmp.2,file = paste0(my.filename[i],".txt"),quote = F,col.names = T,row.names = F,sep="\t")
# }
######################
# real annalysis part#
######################

my.filelist<-list.files(pattern = ".txt")
my.filename<-unlist(lapply(strsplit(my.filelist,"[.]"),"[",1))
for(i in 1:length(my.filelist)){
  assign(my.filename[i],read.table(my.filelist[i],header = T))
  print(paste("finish",my.filelist[i]))
}
# separate condtion and group
my.sample.condition<-read.table("../sample_condition.txt",header = T)
my.sample.condition.2<-my.sample.condition[match(gsub("MAA","",my.filename),my.sample.condition$NAME),]
my.group<-my.sample.condition.2$Disease
my.methylation.cov.list<-list(MAA4130,MAA4308,MAA4382,MAA4414,MAA4494,MAA4617,MAA4660,MAA4788,MAA4807,MAA5190)

my.object<-makeBSseqData(my.methylation.cov.list,my.group)
















