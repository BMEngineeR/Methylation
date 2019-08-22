options(stringsAsFactors = F)
library(Gviz)
library(genomation)
library(dplyr)
library(org.Hs.eg.db)
# change working directory
setwd("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Change_chr/")
# input data
my.data<-read.csv("DM_SNP_all.csv",header = T)
# add filter q-value < 0.05,distance=[-1000,+1000]
my.data.filtered<-filter(my.data,qvalue<=0.05)
my.data.filtered<-my.data.filtered[abs(my.data.filtered$dist.to.feature)<=1000,]

# extract information for creating ned file from whole data. 
my.path.genes<-read.csv("autophagy list.csv")
my.gene.list<-unique(my.path.genes$Approved.symbol) 

my.gene.data<-c()
for (i in 1:length(my.gene.list)){
  gene.tmp<-paste0("^",my.gene.list[i],"$")
  data.tmp<-my.data.filtered[grep(gene.tmp,my.data.filtered$Gene),]
  my.gene.data<-rbind.data.frame(my.gene.data,data.tmp)
  print(i)
}
my.significant.gene<-unique(my.gene.data$Gene)
my.bed<-my.gene.data[,c(1,2,3,7,9)]
# rename colname and rowname
colnames(my.bed)[4]<-"score"
rownames(my.bed)<-1:nrow(my.bed)
head(my.bed)
# create plot function
plot.track<-function(number=1){
  # gene location in bed file
  my.OneGene.bed<-my.bed[grep(paste0("^",my.significant.gene[number],"$"),my.gene.data$Gene),]
  chr.tmp<-unique(my.OneGene.bed$chr)
  # import track information
  itrack<-IdeogramTrack(genome = "hg38",chromosome = unique(my.OneGene.bed$chr))
  # create methylation object for Gviz 
  gtrack<-GenomeAxisTrack() #gene range 
  # gene track
  ucscGenes <- UcscTrack(genome="hg38", table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                         chromosome=chr.tmp, rstarts = "exonStarts", rends = "exonEnds",
                         gene = "name", symbol = 'name', transcript = "name",
                         strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
  z <- ranges(ucscGenes)
  mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
  ucscGenes2 <- ucscGenes
  ranges(ucscGenes2) <- z
  my.ucsc.tmp<-as(ucscGenes2,"data.frame")
  my.gene.track.df<-my.ucsc.tmp[grep(my.significant.gene[number],my.ucsc.tmp$symbol),]
  #my.gene.track.df<-my.gene.track.df[grep("^XM",my.gene.track.df$gene),]
  atrack.tmp<-as(my.OneGene.bed[,c(1:3)],"GRanges")
  atrack<-AnnotationTrack(atrack.tmp,stacking = "dense",fill="white",shape="box")
  # plotTracks(atrack,from = min(my.OneGene.bed$start)-5000,
  #            to=max(my.OneGene.bed$start)+2000)
  # 
  # for(i in 1:nrow(my.gene.track.df)){
  #   if (my.gene.track.df$X.start[i]>my.gene.track.df$X.start[i+1]){
  #     my.gene.track.df<-my.gene.track.df[1:i,]
  #     break()
  #   }
  # }
  my.gene.track.GRange<-as(my.gene.track.df,"GRanges")
  grtrack<-GeneRegionTrack(my.gene.track.GRange,genome = "hg38",
                           chromosome = chr.tmp,
                           name = unique(my.gene.track.df$symbol))
  dtrack<-DataTrack(data=my.OneGene.bed$score,start = my.OneGene.bed$start,end = my.OneGene.bed$end,
                    genome = "hg38", chromosome = chr.tmp,type="histogram",
                    name = "Alzheimer - Normal")
  my.direction<-unique(as.character(my.gene.track.df$X.strand))
  if(my.direction=="-"){
    plotTracks(list(itrack,gtrack,grtrack,atrack,dtrack),from = min(my.OneGene.bed$start)-5000,
               to=max(my.OneGene.bed$start)+2000,type="h",col.main = "black",cex.main = 1.5,
               main = paste(my.significant.gene[number],"TSS visualization"),
               sizes = c(1,1,1,2,2))
  }else{
    plotTracks(list(itrack,gtrack,grtrack,atrack,dtrack),from = min(my.OneGene.bed$start)-2000,
               to=max(my.OneGene.bed$start)+5000,type="mountain",
               type="h",col.main = "black",cex.main = 1.5,
               main = paste(my.significant.gene[number],"TSS visualization"),
               sizes = c(1,1,1,2,2)
               )
  }
}

my.significant.gene
# change output directory. 
# setwd("/BigData/analysis_work/Amer/Amer_run/Sorted_BAM/Change_chr/Gene_plot/")
# for(i in 1:length(my.significant.gene)){
#   tiff(filename = paste0(my.significant.gene[i],"_plot.tiff"))
#   plot.track(i)
#   dev.off()
# }

plot.track(number = 15)






