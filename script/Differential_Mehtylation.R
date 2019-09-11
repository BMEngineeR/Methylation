library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
options(stringsAsFactors = F)
library(Gviz)
library(genomation)
library(dplyr)
library(org.Hs.eg.db)
library(gridExtra)
library(ggthemes)
library(conflicted)
# # preprocessing###########################################################
# suppressMessages(library(methylKit))
# setwd("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/New_analysis_Bowtie2/Sorted_Bam/")
# # read in condition file # commentable
# my.condition<-read.table("sample_condition.txt",header = T)
# my.disease.factor<-factor(my.condition$Disease,labels = c(0,1),levels  = c("Normal","Al"))
# my.sample.id <- as.list(paste0("MAA",my.condition$NAME))
# my.treatment <-  as.numeric(as.character(my.disease.factor))
# # begin with bam 
# # after running this, please comment 

# my.bam.list<- as.list(my.condition$bam)
# my.condition$bam<-paste0("MAA",my.condition$NAME,".sort.bam")
# my.meth<-processBismarkAln(location = my.bam.list,
#                               sample.id = my.sample.id,assembly = "GRCh38",read.context = "CpG",
#                               treatment = my.treatment,
#                               save.folder = getwd())
# remove alternative chromosome,
# after running this, please comment
# my.filelist<-list.files(pattern = ".CpG.txt")
# my.header<-unlist(lapply(strsplit(my.filelist,"_"),'[[',1))
# for (i in 1:length(my.filelist)){
#   x<-read.table(my.filelist[i],header = T)
#   x$chrBase<-paste0("chr",x$chrBase)
#   x$chr<-paste0("chr",x$chr)
#   xx<-x[grep("^chr[1-9,X,Y,MT]",x$chr),]
#   write.table(xx,
#               file = paste0(my.header[i],"_CpGs.chr.txt"),
#               quote = F,
#               row.names = F,
#               sep = "\t")
#   print(i)
# }
# # read in data
#####################################################################
suppressMessages(library(methylKit))
setwd("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/New_analysis_Bowtie2/Sorted_Bam/")
# read in condition file # commentable
my.condition<-read.table("sample_condition.txt",header = T)
# make name consistency # commentable
rownames(my.condition)<-paste0("MAA",my.condition$NAME,"_CpGs.chr.txt")
# make order consistency # commentable
# my.condition<-my.condition[match(list.files(pattern = "_CpGs.chr.txt"),row.names(my.condition)),]
my.disease.factor<-factor(my.condition$Disease,labels = c(0,1),levels  = c("Normal","Al"))
# creating Meth object.
my.filelist<-as.list(rownames(my.condition))
my.sample.id <- as.list(paste0("MAA",my.condition$NAME))
my.treatment <-  as.numeric(as.character(my.disease.factor))
my.meth<-methRead(my.filelist,
                  sample.id=my.sample.id,
                  assembly="GRCh38",
                  treatment=my.treatment,
                  context="CpG")
save(my.meth,file = "my.meth.rds")
load("my.meth.rds")
# visually showing the methylation status
# par(mfrow=c(2,2))
# for(i in 1:4){
#   getMethylationStats(my.meth[[i]],plot = T,both.strands = F)
# }
# getMethylationStats(my.meth[[1]],plot = T,both.strands = F)
# getCoverageStats(my.meth[[1]],plot = T,both.strands =F)
# # Filtering samples based on the read coverage

my.filterd.meth<-filterByCoverage(my.meth,lo.count = 10,lo.perc = NULL,
                                  hi.count = NULL, hi.perc = 99.9)
# merger all samples
my.meth.merged <- unite(my.filterd.meth,destrand = FALSE)

# look at sample correlation
# pdf("correlationPlot.pdf")
# getCorrelation(my.meth.merged,plot = TRUE)
# dev.off()
# 
# # clutering samples and plot PCA
head(my.meth.merged)
clusterSamples(my.meth.merged,dist = "correlation", method = "ward.D", plot = TRUE)

PCASamples(my.meth.merged)

# Find differentially methylated bases or regions
MyDiff <- calculateDiffMeth(my.meth.merged)
# select significantly differential region. 
MyDiff.significant<-getMethylDiff(MyDiff,difference=5,qvalue=0.05)
# save file for later analysis

# save differential Methylation file for checkpoint. 
save(MyDiff,file = "MyDiff.rds")
save(MyDiff.significant,file = "MyDiff_significant.rds")
###################################################
#loading data#####################################
load("MyDiff.rds")
load("MyDiff_significant.rds")
# annotation
library(genomation)
gene.obj=readTranscriptFeatures("GENE.location.Info.txt")
diff_gene<-annotateWithGeneParts(as(MyDiff,"GRanges"),gene.obj)

my.Diff.anno<-diff_gene@dist.to.TSS
my.Diff.meth<-as.data.frame(MyDiff@.Data)
colnames(my.Diff.meth)<-c("chr","start","end","strand","pvalue",'qvalue',"difference")
my.Diff.meth.clean<-my.Diff.meth[grep("^chr[1-9,X,Y,MT]",my.Diff.meth$chr),]

my.Diff.final<- cbind(my.Diff.meth.clean,my.Diff.anno[,c(2,3)])
gene.refseq.name<-unlist(sapply(strsplit(my.Diff.final$feature.name,"[.]"),'[[',1))
gene.Symbol<-AnnotationDbi::select(org.Hs.eg.db,keys = gene.refseq.name,
                                   keytype = "REFSEQ",
                                   columns = "SYMBOL")
x.new<-cbind.data.frame(my.Diff.final,feature.name=gene.refseq.name,Gene=gene.Symbol$SYMBOL)
# generatecsv
write.csv(x.new,file = "DM_SNP_all.csv",quote = F,row.names = F)


# group 1(treatment) - group 0(control) : if treatment > control, then difference in  my.Diff.final should be postive number. 
#################################################################
#####quality control and general situation#######################
#################################################################
diff_gene<-annotateWithGeneParts(as(MyDiff,"GRanges"),gene.obj)
head(getAssociationWithTSS(diff_gene))
# show basic statisitc
getTargetAnnotationStats(diff_gene,percentage=TRUE,precedence=TRUE)
getTargetAnnotationStats(diff_gene,percentage = T,precedence = T)
# plot each region percentage
plotTargetAnnotation(diff_gene,precedence=TRUE,
                     main="differential methylation annotation")


hist(my.Diff.final$dist.to.feature,breaks = 1000,col = "skyblue",main = "Distance to TSS",xlim = c(-100000,100000),xlab = "distance")

