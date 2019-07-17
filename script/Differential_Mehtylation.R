suppressMessages(library(methylKit))
setwd("/BigData/analysis_work/Amer/Amer_run/Sorted_BAM/Change_chr/")
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
# my.meth<-processBismarkAln(location = my.filelist,
#                               sample.id = my.sample.id,assembly = "GRCh38",read.context = "CpG",
#                               treatment = my.treatment,
#                               save.folder = getwd())
my.meth<-methRead(my.filelist,
                  sample.id=my.sample.id,
                  assembly="GRCh38",
                  treatment=my.treatment,
                  context="CpG"
)
# visually showing the methylation status
par(mfrow=c(2,2))
for(i in 1:4){
  getMethylationStats(my.meth[[i]],plot = T,both.strands = F)
}
getMethylationStats(my.meth[[1]],plot = T,both.strands = F)
getCoverageStats(my.meth[[1]],plot = T,both.strands =F)
# Filtering samples based on the read coverage

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
MyDiff.all<-getMethylDiff(MyDiff,difference=10,qvalue=0.05)
write.table(MyDiff.all,file="my.diff.meth.txt",row.names = F,quote = F)
# save differential Methylation file for checkpoint. 
save(MyDiff,file = "MyDiff.rds")
# visualize distribution of hypo-meth/hypo-methylated base
Diff.Meth.report<-diffMethPerChr(MyDiff.all,plot=F,cutoff=0.05,meth.cutoff = 10)
write.csv(Diff.Meth.report,file = "Diff.Meth.Report.csv",quote = F,row.names = F )
# annotation
library(genomation)
gene.obj=readTranscriptFeatures("GENE.location.Info.txt")
diff_gene<-annotateWithGeneParts(as(MyDiff,"GRanges"),gene.obj)
# cpg.obj<-readFeatureFlank("hgTables.CpGi.txt",feature.flank.name=c("CpGi","shores"))
# diffCpGann<-annotateWithFeatureFlank(as(MyDiff.all,"GRanges"),
#                                     cpg.obj$CpGi,cpg.obj$shores,
#                                     feature.name="CpGi",flank.name="shores")
# create differential matrix(SNPs)
# for(i in 1:dim(my.Diff.anno)[1]){
#   tmp<-my.Diff.anno$target.row[i+1]-my.Diff.anno$target.row[i]
#   if(tmp==1){print(i)} else{break()}
# }
my.Diff.anno<-diff_gene@dist.to.TSS
my.Diff.meth<-as.data.frame(MyDiff@.Data)
colnames(my.Diff.meth)<-c("chr","start","end","strand","pvalue",'qvalue',"difference")
my.Diff.meth.clean<-my.Diff.meth[grep("^chr[1-9,X,Y,MT]",my.Diff.meth$chr),]

my.Diff.final<- cbind(my.Diff.meth.clean,my.Diff.anno[,c(2,3)])
# group 1(treatment) - group 0(control) : if treatment > control, then difference in  my.Diff.final should be postive number. 

#write.table(diff_gene@dist.to.TSS,file = "DM.txta")
head(getAssociationWithTSS(diff_gene))
getTargetAnnotationStats(diff_gene,percentage = T,precedence = T)
plotTargetAnnotation(diff_gene,precedence=TRUE,
                     main="differential methylation annotation")
# debug for what hell happen?? why the row has different number? 
# because they introduced the alternative chromesome, no annotation file for that. 
# for(i in 19917:21176){
#   tmp_d1<-my.Diff.meth$start[i+1]-my.Diff.meth$start[i]
#   tmp_d2<-my.Diff.anno$dist.to.feature[i+1]- my.Diff.anno$dist.to.feature[i]
#   tmp_d3<-my.Diff.anno$target.row[i+1]- my.Diff.anno$target.row[i]
#   if (abs(tmp_d3)==1){
#     print(paste("ok for",i))
#   } else{
#     print(paste("warning for",i))
#     break()
#     }
# }
## pathway enrichment.
setwd("/BigData/analysis_work/Amer/Amer_run/Sorted_BAM/Change_chr/report/")
# my.Diff.final<-read.csv("DM_SNP.csv",header = T,stringsAsFactors = F)
hist(my.Diff.final$dist.to.feature,breaks = 1000,col = "skyblue",main = "Distance to TSS",xlim = c(-100000,100000),xlab = "distance")
# my.promoter.table<-my.Diff.final[which(abs(my.Diff.final$dist.to.feature)<=1000),]
# treat.promoter.table<-my.promoter.table[which(my.promoter.table$difference>0),]
# treat.promoter.gene<-as.character(treat.promoter.table$feature.name)
