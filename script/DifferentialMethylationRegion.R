if (!require(DSS)){
  BiocManager::install("DSS")
}
require(bsseq)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(genomation)
options(stringsAsFactors = F)
setwd("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/coverage_meth")
# prepressing file, they contain alternative chromosome
my.filelist<-list.files(pattern = "cov.gz")
my.filename<-unlist(lapply(strsplit(my.filelist,"_"),"[",1))
for (i in (1:length(my.filelist))){
  #read in file
  tmp<-read.table(my.filelist[i])
  # filter out alternative chromosome
  tmp.1<-tmp[grep("^[1-9,X,Y,M]",tmp$V1),]
  # add chrX instead of X
  tmp.1[,1]<-paste0("chr",tmp.1$V1)
  tmp.2<-cbind.data.frame(chr=tmp.1$V1,pos=tmp.1$V2,N=c(tmp.1$V5+tmp.1$V6),X=tmp.1$V5)
  write.table(tmp.2,file = paste0(my.filename[i],".txt"),quote = F,col.names = T,row.names = F,sep="\t")
}
######################
# real annalysis part#
######################

my.filelist<-list.files(pattern = "[0-9].txt")
my.filename<-unlist(lapply(strsplit(my.filelist,"[.]"),"[",1))
for(i in 1:length(my.filelist)){
  assign(my.filename[i],read.table(my.filelist[i],header = T))
  print(paste("finish",my.filelist[i]))
}
# separate condtion and group
my.sample.condition<-read.table("../sample_condition.txt",header = T)
my.sample.condition$Group<-paste0(my.sample.condition$Disease,".",my.sample.condition$NAME)
my.sample.condition.2<-my.sample.condition[match(gsub("MAA","",my.filename),my.sample.condition$NAME),]
my.group<-my.sample.condition.2$Group
my.methylation.cov.list<-list(MAA4130,MAA4308,MAA4382,MAA4414,MAA4494,MAA4617,MAA4660,MAA4788,MAA4807,MAA5190)
# checkpoint
# save(my.methylation.cov.list,file = "my.methylation.cov.list")
# load("my.methylation.cov.list")
my.object<-makeBSseqData(my.methylation.cov.list,my.group)
# save(my.object,file="my.object")
 load("./my.object")
# RRBS does not recommend smoothing process
my.normal<-my.group[grep("^Normal",my.group)]
my.alzheimer<-my.group[grep("^Al",my.group)]
DMLTest<- DMLtest(my.object,group1 = my.normal,group2 = my.alzheimer)
#############################
#############checkpoint######
#############################
dmls <-callDML(DMLTest,p.threshold = 0.05)
save(dmls,file = 'dmls')
Sys.time()
dmrs <- callDMR(DMLTest,p.threshold = 0.05)

###########################
####  check point##########
###########################
save(dmrs,file = "dmrs")##
load("dmrs")            ##
###########################
head(dmrs)
# annotate gene 
gene.obj=readTranscriptFeatures("GENE.location.Info.txt")
gene.annotate<-annotateWithGeneParts(as(dmrs,"GRanges"),gene.obj)
head(gene.annotate@dist.to.TSS)
# add annotated gene to new.DMR.v1
new.DMR.v1<-cbind.data.frame(dmrs,gene.annotate@dist.to.TSS[,c(2,3,4)],gene.annotate@members)
# extract ref. name from new.DMR.v1
gene.REFSEQ<-gsub(".[0-9]+$","",new.DMR.v1$feature.name)
# switch name from refseq to gene symbol
gene.SYMBOL<-AnnotationDbi::select(org.Hs.eg.db,keys = gene.REFSEQ,
                                   keytype = "REFSEQ",
                                   columns = "SYMBOL")
# insery gene symbol to new DMR
new.DMR.v2<-cbind.data.frame(new.DMR.v1[,1:11],gene.symbol=gene.SYMBOL$SYMBOL,new.DMR.v1[,12:15])
write.csv(new.DMR.v2,file="DMR.csv",quote = F,row.names = F)









