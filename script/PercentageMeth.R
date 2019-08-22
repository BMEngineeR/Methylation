options(stringsAsFactors = F)
library(Gviz)
library(genomation)
library(dplyr)
library(org.Hs.eg.db)
library(gridExtra)
library(ggthemes)
suppressMessages(library(methylKit))
setwd("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Change_chr/")
# read in condition file # commentable
my.condition<-read.table("sample_condition.txt",header = T)
# make name consistency # commentable
rownames(my.condition)<-paste0("MAA",my.condition$NAME,"_CpGs.chr.txt")
# make order consistency # commentable
my.condition<-my.condition[match(list.files(pattern = "_CpGs.chr.txt"),row.names(my.condition)),]
my.disease.factor<-factor(my.condition$Disease,labels = c(0,1),levels  = c("Normal","Al"))
# creating Meth object.
my.filelist<-as.list(list.files(pattern = "CpGs.chr.txt"))
my.sample.id <- as.list(paste0("MAA",my.condition$NAME))
my.treatment <-  as.numeric(as.character(my.disease.factor))
my.meth<-methRead(my.filelist,
                  sample.id=my.sample.id,
                  assembly="GRCh38",
                  treatment=my.treatment,
                  context="CpG"
)
my.filterd.meth<-filterByCoverage(my.meth,lo.count = 10,lo.perc = NULL,
                                  hi.count = NULL, hi.perc = 99.9)
my.meth.normalized<-normalizeCoverage(my.filterd.meth)
my.meth.merged <- unite(my.meth.normalized,destrand = F)
# MyDiff <- calculateDiffMeth(my.meth.merged)
# # select significantly differential region. 
# MyDiff.all<-getMethylDiff(MyDiff,difference=10,qvalue=0.05)
my.meth.merged<-as.data.frame(my.meth.merged)
# read in all data
my.data<-read.csv("DM_SNP_all.csv",header = T)
# add filter q-value < 0.05,distance=[-1000,+1000]
my.data.filtered<-filter(my.data,qvalue<=0.05)
my.data.filtered<-my.data.filtered[abs(my.data.filtered$dist.to.feature)<=1000,]
# before gene my.significant.gene please run GenomeVisulization.R
my.significant.gene

plot.gene<-function(gene.name="ATG5",...){
  my.target.gene<-my.data.filtered[grep(gene.name,my.data.filtered$Gene),]
  my.chr<-unique(my.target.gene$chr)
  my.start<-my.target.gene$start
  
  
  my.merged.FilterByChr<-my.meth.merged[grep(paste0("^",my.chr,"$"),my.meth.merged$chr),]
  my.merged.FilterByStart<-sapply(my.start,function(x) my.merged.FilterByChr[grep(paste0("^",x,"$"),my.merged.FilterByChr$start),])
  my.merged.FilterByStart<-t(as.data.frame(my.merged.FilterByStart))
  my.condition.table<-cbind.data.frame(ID=my.meth.merged@sample.ids,Condition=my.meth.merged@treatment)
  my.condition.table.sort<-my.condition.table[order(my.condition.table$Condition),]
  match.index<-match(my.condition.table.sort$ID,my.condition.table$ID)
  # match.index.new<-c(1,2,5,6,7,8,9,10,17,18,19,23,24,25,32,33,34,11,12,13,14,15,16,20,21,22,26,27,28,29,30,31)
  # my.plot.table<-my.merged.FilterByStart[,match.index.new]
  normal.mean.coverage <-rowMeans(matrix(as.numeric(my.merged.FilterByStart[,3*match.index[1:5]+2]),ncol = 5))
  normal.mean.Cs<-rowMeans(matrix(as.numeric(my.merged.FilterByStart[,3*match.index[1:5]+3]),ncol = 5))
  al.mean.coverage <-rowMeans(matrix(as.numeric(my.merged.FilterByStart[,3*match.index[6:10]+2]),ncol = 5))
  al.mean.Cs<-rowMeans(matrix(as.numeric(my.merged.FilterByStart[,3*match.index[6:10]+3]),ncol = 5))
  normal.percent.Cs<-(normal.mean.Cs/normal.mean.coverage)*100
  al.percent.Cs<-(al.mean.Cs/al.mean.coverage)*100
  
  my.plot.df<-cbind.data.frame(chr=as.array(as.character(my.merged.FilterByStart[,1])),
                               start=as.array(as.numeric(my.merged.FilterByStart[,2])),
                               al.mean.Cs=al.mean.Cs,
                               al.mean.Ts=al.mean.coverage-al.mean.Cs,
                               al.mean.coverage=al.mean.coverage,
                               al.per=al.percent.Cs,
                               normal.mean.Cs=normal.mean.Cs,
                               normal.mean.Ts=normal.mean.coverage- normal.mean.Cs,
                               normal.mean.coverage=normal.mean.coverage,
                               normal.per=normal.percent.Cs)
  my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
  #my.plot.df$y<-rep(100,nrow(my.plot.df))
  
  p.normal <- ggplot(my.plot.df,aes(x=start,y=al.per))+
    geom_segment(aes(x=start,y=0,xend=start,yend=100),color="#20639B",size=1,data=my.plot.df)+
    geom_segment(aes(x=start,y=0,xend=start,yend=normal.per),color="#C70039",size=1,data=my.plot.df)+
    theme(panel.background = element_rect(fill="white", colour="white", size=0.1, 
                                          linetype="blank", color="white"),
          axis.text.x = element_blank())+
    geom_hline(yintercept=c(0,50,100), linetype="solid", color="black", size=0.5)+
    labs(title=paste0(gene.name," Gene Methylation"),
         x ="", y = "Normal")
  
  
  p.Al <- ggplot(my.plot.df,aes(x=start,y=al.per))+
    geom_segment(aes(x=start,y=0,xend=start,yend=100),color="#20639B",size=1,data=my.plot.df)+
    geom_segment(aes(x=start,y=0,xend=start,yend=al.per),color="#C70039",size=1,data=my.plot.df)+
    theme(panel.background = element_rect(fill="white", colour="white", size=0.1, 
                                          linetype="blank", color="white"))+
    geom_hline(yintercept=c(0,50,100), linetype="solid", color="black", size=0.5)+
    labs(title="",
         x ="", y = "Alzheimer")
  
  grid.arrange(p.normal,p.Al)
}
my.significant.gene
plot.gene(gene.name = "ATG16L2")













## coverage check
my.sample.name<-my.meth.merged@sample.ids
my.readnumber<-c(33171523,43165280,29397200,46211060,39995359,
                      25921837,30249793,31428974,35737181,31702113)
my.coverage<-150*my.readnumber/3.2/10^9
my.report<-rbind(as.character(my.readnumber),as.character(my.coverage))
colnames(my.report)<-my.sample.name
rownames(my.report)<-c("reads number","coverage")
my.report


