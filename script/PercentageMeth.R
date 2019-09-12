options(stringsAsFactors = F)
library(Gviz)
library(genomation)
library(dplyr)
library(org.Hs.eg.db)
library(gridExtra)
library(ggthemes)
library(ggplot2)
suppressMessages(library(methylKit))
setwd("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/New_analysis_Bowtie2/Sorted_Bam/")
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
my.filterd.meth<-filterByCoverage(my.meth,lo.count = NULL,lo.perc = NULL,
                                  hi.count = NULL, hi.perc = 99.9)
my.meth.normalized<-normalizeCoverage(my.filterd.meth)
my.meth.merged <- unite(my.meth.normalized,destrand = F)
# MyDiff <- calculateDiffMeth(my.meth.merged)
# # select significantly differential region. W
# MyDiff.all<-getMethylDiff(MyDiff,difference=10,qvalue=0.05)
my.meth.merged<-as.data.frame(my.meth.merged)
# read in all data
my.data<-read.csv("DM_SNP_all.csv",header = T)
# add filter q-value < 0.05,distance=[-1000,+1000]
my.data.filtered<-dplyr::filter(my.data,qvalue<=0.05)
my.data.filtered<-my.data.filtered[abs(my.data.filtered$dist.to.feature)<=1000,]
# before gene my.significant.gene please run GenomeVisulization.R
my.significant.gene

plot.gene<-function(gene.name="ATG5",...){
  my.target.gene<-my.data.filtered[grep(paste0("^",gene.name,"$"),my.data.filtered$Gene),]
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

gene.region.list<-read.table("GeneRegion.txt")
colnames(gene.region.list)<- c("name","region")

gene.region=gene.region.list

  
plot.region<-function(Gene.Region=gene.region,index=2,...){
  # split gene region to chrome, start, end
  tmp.gene.info<-Gene.Region[index,]
  # gene name
  tmp.gene.name<-as.character(tmp.gene.info[1])
  tmp.gene.region<-as.character(tmp.gene.info[2])
  # chromosome
  tmp.gene.chr<-strsplit(tmp.gene.region,":")[[1]][1]
  tmp.gene.StartToEnd<-strsplit(tmp.gene.region,":")[[1]][2]
  # start
  tmp.gene.start<-strsplit(tmp.gene.StartToEnd,"-")[[1]][1]
  # end
  tmp.gene.end<-strsplit(tmp.gene.StartToEnd,"-")[[1]][2]
  
  my.merged.FilterByChr<-my.meth.merged[grep(paste0("^",tmp.gene.chr,"$"),my.meth.merged$chr),]
  tmp.gene.start.iter<-as.numeric(tmp.gene.start)
  while(!tmp.gene.start.iter %in% my.merged.FilterByChr$start){
    tmp.gene.start.iter<-tmp.gene.start.iter-1
    # print(tmp.gene.start.iter)
  }
  
  tmp.gene.start.index<-grep(tmp.gene.start.iter,my.merged.FilterByChr$start)
  tmp.gene.end.iter<-as.numeric(tmp.gene.end)
  while(!tmp.gene.end.iter %in% my.merged.FilterByChr$start){
    tmp.gene.end.iter<-tmp.gene.end.iter+1
     # print(tmp.gene.end.iter)
  } 
  # get location of region end in my.meth.merged 
  tmp.gene.end.index<-grep(tmp.gene.end.iter,my.merged.FilterByChr$start)
  # region index is aim location 
  region.index<-seq(tmp.gene.start.index,tmp.gene.end.index)
  my.merged.FilterByStart<-my.merged.FilterByChr[region.index,]
  
  my.condition.table<-cbind.data.frame(ID=my.meth.merged@sample.ids,Condition=my.meth.merged@treatment)
  my.condition.table.sort<-my.condition.table[order(my.condition.table$Condition),]
  match.index<-match(my.condition.table.sort$ID,my.condition.table$ID)
  # match.index.new<-c(1,2,5,6,7,8,9,10,17,18,19,23,24,25,32,33,34,11,12,13,14,15,16,20,21,22,26,27,28,29,30,31)
  # my.plot.table<-my.merged.FilterByStart[,match.index.new]
  normal.mean.coverage <-rowMeans(getData(my.merged.FilterByStart)[,3*match.index[1:5]+2])
  normal.mean.Cs<-rowMeans(getData(my.merged.FilterByStart)[,3*match.index[1:5]+3])
  al.mean.coverage <-rowMeans(getData(my.merged.FilterByStart)[,3*match.index[6:10]+2])
  al.mean.Cs<-rowMeans(getData(my.merged.FilterByStart)[,3*match.index[6:10]+3])
  normal.percent.Cs<-(normal.mean.Cs/normal.mean.coverage)*100
  al.percent.Cs<-(al.mean.Cs/al.mean.coverage)*100
  print("ploting...")
  my.plot.df<-cbind.data.frame(chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                               start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
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
  
  p.normal <- ggplot(my.plot.df,aes(x=start,y=normal.per))+
    geom_segment(aes(x=start,y=0,xend=start,yend=100),color="#20639B",size=1,data=my.plot.df)+
    geom_segment(aes(x=start,y=0,xend=start,yend=normal.per),color="#C70039",size=1,data=my.plot.df)+
    theme(panel.background = element_rect(fill="white", colour="white", size=0.1, 
                                          linetype="blank", color="white"),
          plot.title = element_text(size=11),
          axis.text.x = element_blank())+
    geom_hline(yintercept=c(0,50,100), linetype="solid", color="black", size=0.5)+
    labs(title=paste0(tmp.gene.name," Gene Methylation",
                      ": (",tmp.gene.chr," : ",tmp.gene.start,"-",tmp.gene.end," )"),
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
plot.region(Gene.Region = gene.region,index = 2)

plot.region.separate<-function(Gene.Region=gene.region,index=2,...){
  # split gene region to chrome, start, end
  tmp.gene.info<-Gene.Region[index,]
  # gene name
  tmp.gene.name<-as.character(tmp.gene.info[1])
  tmp.gene.region<-as.character(tmp.gene.info[2])
  # chromosome
  tmp.gene.chr<-strsplit(tmp.gene.region,":")[[1]][1]
  tmp.gene.StartToEnd<-strsplit(tmp.gene.region,":")[[1]][2]
  # start
  tmp.gene.start<-strsplit(tmp.gene.StartToEnd,"-")[[1]][1]
  # end
  tmp.gene.end<-strsplit(tmp.gene.StartToEnd,"-")[[1]][2]
  
  my.merged.FilterByChr<-my.meth.merged[grep(paste0("^",tmp.gene.chr,"$"),my.meth.merged$chr),]
  tmp.gene.start.iter<-as.numeric(tmp.gene.start)
  while(!tmp.gene.start.iter %in% my.merged.FilterByChr$start){
    tmp.gene.start.iter<-tmp.gene.start.iter-1
    # print(tmp.gene.start.iter)
  }
  
  tmp.gene.start.index<-grep(tmp.gene.start.iter,my.merged.FilterByChr$start)
  tmp.gene.end.iter<-as.numeric(tmp.gene.end)
  while(!tmp.gene.end.iter %in% my.merged.FilterByChr$start){
    tmp.gene.end.iter<-tmp.gene.end.iter+1
    # print(tmp.gene.end.iter)
  } 
  # get location of region end in my.meth.merged 
  tmp.gene.end.index<-grep(tmp.gene.end.iter,my.merged.FilterByChr$start)
  # region index is aim location 
  region.index<-seq(tmp.gene.start.index,tmp.gene.end.index)
  my.merged.FilterByStart<-my.merged.FilterByChr[region.index,]
  
  my.condition.table<-cbind.data.frame(ID=my.meth.merged@sample.ids,Condition=my.meth.merged@treatment)
  my.condition.table.sort<-my.condition.table[order(my.condition.table$Condition),]
  match.index<-match(my.condition.table.sort$ID,my.condition.table$ID)
  # match.index.new<-c(1,2,5,6,7,8,9,10,17,18,19,23,24,25,32,33,34,11,12,13,14,15,16,20,21,22,26,27,28,29,30,31)
  # my.plot.table<-my.merged.FilterByStart[,match.index.new]
  normal.sample.names<-my.merged.FilterByStart@sample.ids[match.index[1:5]]
  al.sample.names<-my.merged.FilterByStart@sample.ids[match.index[6:10]]
  normal.all.coverage <-getData(my.merged.FilterByStart)[,3*match.index[1:5]+2]
  normal.all.Cs<-getData(my.merged.FilterByStart)[,3*match.index[1:5]+3]
  al.all.coverage <-getData(my.merged.FilterByStart)[,3*match.index[6:10]+2]
  al.all.Cs<-getData(my.merged.FilterByStart)[,3*match.index[6:10]+3]
  normal.percent.Cs<-(normal.all.Cs/normal.all.coverage)*100
  al.percent.Cs<-(al.all.Cs/al.all.coverage)*100
  # load function
  plot.normal<-function(plot.df=NULL,name=NULL){
    p.normal <- ggplot(plot.df,aes(x=start,y=normal.per))+
      geom_segment(aes(x=start,y=0,xend=start,yend=100),color="#20639B",size=1,data=my.plot.df)+
      geom_segment(aes(x=start,y=0,xend=start,yend=normal.per),color="#C70039",size=1,data=my.plot.df)+
      theme(panel.background = element_rect(fill="white", colour="white", size=0.1, 
                                            linetype="blank", color="white"),
            plot.title = element_text(size=11),axis.title.y = element_text(size = 10),
            axis.text.x = element_blank())+
      geom_hline(yintercept=c(0,50,100), linetype="solid", color="black", size=0.5)+
      labs(title=NULL,
           x ="", y = paste0("Normal:",name))
    return(p.normal)
  }
  # sample1 normal1
    sample1.name<-normal.sample.names[1]
    my.plot.df<-cbind.data.frame(name=rep(sample1.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 normal.per=(normal.all.Cs[,1]/normal.all.coverage[,1])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Nor.sample1.plot<-   ggplot(my.plot.df,aes(x=start,y=normal.per))+
      geom_segment(aes(x=start,y=0,xend=start,yend=100),color="#20639B",size=1,data=my.plot.df)+
      geom_segment(aes(x=start,y=0,xend=start,yend=normal.per),color="#C70039",size=1,data=my.plot.df)+
      theme(panel.background = element_rect(fill="white", colour="white", size=0.1, 
                                            linetype="blank", color="white"),
            plot.title = element_text(size=11),axis.title.y = element_text(size = 10),
            axis.text.x = element_blank())+
      geom_hline(yintercept=c(0,50,100), linetype="solid", color="black", size=0.5)+
      labs(title=paste0(tmp.gene.name," Gene Methylation",
                        ": (",tmp.gene.chr," : ",tmp.gene.start,"-",tmp.gene.end," )"),
           x ="", y = paste0("Normal:",sample1.name))
  # sample4 normal4
    sample4.name<-normal.sample.names[4]
    my.plot.df<-cbind.data.frame(name=rep(sample4.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 normal.per=(normal.all.Cs[,4]/normal.all.coverage[,4])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Nor.sample4.plot<-plot.normal(plot.df=my.plot.df,name =sample4.name )
  # sample2 normal2
    sample2.name<-normal.sample.names[2]
    my.plot.df<-cbind.data.frame(name=rep(sample2.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 normal.per=(normal.all.Cs[,2]/normal.all.coverage[,2])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Nor.sample2.plot<-plot.normal(plot.df=my.plot.df,name =sample2.name)
    # sample3 normal3 
    sample3.name<-normal.sample.names[3]
    my.plot.df<-cbind.data.frame(name=rep(sample3.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 normal.per=(normal.all.Cs[,3]/normal.all.coverage[,3])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Nor.sample3.plot<-plot.normal(plot.df=my.plot.df,name =sample3.name)
    # sample5 normal5
    sample5.name<-normal.sample.names[5]
    
    my.plot.df<-cbind.data.frame(name=rep(sample5.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 normal.per=(normal.all.Cs[,5]/normal.all.coverage[,5])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Nor.sample5.plot<-plot.normal(plot.df=my.plot.df,name =sample5.name)
# plot alzheimer
    plot.al<-function(plot.df=NULL,name=NULL){
      p.Al <- ggplot(plot.df,aes(x=start,y=al.per))+
        geom_segment(aes(x=start,y=0,xend=start,yend=100),color="#20639B",size=1,data=my.plot.df)+
        geom_segment(aes(x=start,y=0,xend=start,yend=al.per),color="#C70039",size=1,data=my.plot.df)+
        theme(panel.background = element_rect(fill="white", colour="white", size=0.1, 
                                              linetype="blank", color="white"),axis.text.x = element_blank(),axis.title.y = element_text(size = 10))+
        geom_hline(yintercept=c(0,50,100), linetype="solid", color="black", size=0.5)+
        labs(title="",
             x ="", y = paste0("Alzheimer:",name))
      return(p.Al)
    }
    # sample1 Al1
    sample1.name<-al.sample.names[1]
    my.plot.df<-cbind.data.frame(name=rep(sample1.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 al.per=(al.all.Cs[,1]/al.all.coverage[,1])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Al.sample1.plot<-plot.al(plot.df=my.plot.df,name =sample1.name)
    # sample2 Al2
    sample2.name<-al.sample.names[2]
    my.plot.df<-cbind.data.frame(name=rep(sample2.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 al.per=(al.all.Cs[,2]/al.all.coverage[,2])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Al.sample2.plot<-plot.al(plot.df=my.plot.df,name =sample2.name)
    # sample3 Al3
    sample3.name<-al.sample.names[3]
    my.plot.df<-cbind.data.frame(name=rep(sample3.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 al.per=(al.all.Cs[,3]/al.all.coverage[,3])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Al.sample3.plot<-plot.al(plot.df=my.plot.df,name =sample3.name)
    # sample4 Al4
    sample4.name<-al.sample.names[4]
    my.plot.df<-cbind.data.frame(name=rep(sample4.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 al.per=(al.all.Cs[,4]/al.all.coverage[,4])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Al.sample4.plot<-plot.al(plot.df=my.plot.df,name =sample4.name)
    # sample5 Al5
    sample5.name<-al.sample.names[5]
    my.plot.df<-cbind.data.frame(name=rep(sample5.name,nrow(my.merged.FilterByStart)),
                                 chr=as.array(as.character(getData(my.merged.FilterByStart)[,1])),
                                 start=as.array(as.numeric(getData(my.merged.FilterByStart)[,2])),
                                 al.per=(al.all.Cs[,5]/al.all.coverage[,5])*100)
    my.plot.df$x<-my.plot.df$start- min(my.plot.df$start)+1
    Al.sample5.plot<-ggplot(my.plot.df,aes(x=start,y=al.per))+
      geom_segment(aes(x=start,y=0,xend=start,yend=100),color="#20639B",size=1,data=my.plot.df)+
      geom_segment(aes(x=start,y=0,xend=start,yend=al.per),color="#C70039",size=1,data=my.plot.df)+
      theme(panel.background = element_rect(fill="white", colour="white", size=0.1, 
                                            linetype="blank", color="white"),axis.title.y = element_text(size = 10))+
      geom_hline(yintercept=c(0,50,100), linetype="solid", color="black", size=0.5)+
      labs(title="",
           x ="", y = paste0("Alzheimer:",sample5.name))
  
  grid.arrange(Nor.sample1.plot,Nor.sample2.plot,Nor.sample3.plot,Nor.sample4.plot,Nor.sample5.plot,
               Al.sample1.plot,Al.sample2.plot,Al.sample3.plot,Al.sample4.plot,Al.sample5.plot,
               nrow=10)
}
plot.region.separate(Gene.Region = gene.region,index = 3)

for (i in 1:nrow(gene.region)){
  svg(paste0("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/New_analysis_Bowtie2/Sorted_Bam/separate_plot/",gene.region$name[i],"CpGs_Meth.svg"),height = 1200,width = 800)
  plot.region.separate(Gene.Region = gene.region,index =9)
  dev.off()
}

my.significant.gene
plot.gene(gene.name = "ULK2")

my.significant.gene2<-c(my.significant.gene,"NBR1",'DNMT3A','DNMT3B')


for (i in my.significant.gene2){
  svg(paste0("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Change_chr/plot/",i,"CpGs_Meth.svg"),height = 10,width = 10)
  plot.gene(gene.name = i)
  dev.off()
}

plot.gene(gene.name ="ATG5")






## coverage check
my.sample.name<-my.meth.merged@sample.ids
my.readnumber<-c(33171523,43165280,29397200,46211060,39995359,25921837,30249793,31428974,35737181,31702113)
my.readnumber<-2*my.readnumber
my.coverage<-150*my.readnumber/3.2/10^9


my.report<-rbind(as.character(my.readnumber),as.character(my.coverage))
colnames(my.report)<-my.sample.name
rownames(my.report)<-c("reads number","coverage")
my.report


