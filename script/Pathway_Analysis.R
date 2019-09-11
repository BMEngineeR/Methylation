# GO pathway 
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
# setworking directory
setwd("/fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/New_analysis_Bowtie2/Sorted_Bam/")
### read in file.
x.new<-read.csv("DM_SNP_all.csv",header = T,stringsAsFactors = F)
# add filter qvalue<0.05
x.new<-x.new[x.new$qvalue<=0.05,]
# select feature.
# up regulate
# add filter, differen > 10, alzheimer significant 
my.up.coordinate<-intersect(which(x.new$difference>= 10),which(abs(x.new$dist.to.feature)<=1000))
my.gene.up<-unique(x.new$Gene[my.up.coordinate])
my.Go.CC.up<- enrichGO(gene=my.gene.up,OrgDb = org.Hs.eg.db,
                    ont="CC",keyType = "SYMBOL",pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

my.Go.BP.up<- enrichGO(gene=my.gene.up,OrgDb = org.Hs.eg.db,
                    ont="BP",keyType = "SYMBOL",pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

my.Go.MF.up<- enrichGO(gene=my.gene.up,OrgDb = org.Hs.eg.db,
                    ont="MF",keyType = "SYMBOL",pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
my.entrize.id.up<-AnnotationDbi::select(org.Hs.eg.db,keys = my.gene.up,
                                        keytype = "SYMBOL",
                                        columns = "ENTREZID")
my.KEGG.up<-enrichKEGG(gene=my.entrize.id.up$ENTREZID,organism = "hsa",
                       keyType = "kegg",pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
write.table(my.Go.CC.up@result,file = "GO_enrichment_CC_AL.csv",quote = F,row.names = F,sep = ";")
write.table(my.Go.MF.up@result,file = "GO_enrichment_MF_AL.csv",quote = F,row.names = F,sep=";")
write.table(my.Go.BP.up@result,file = "GO_enrichment_BP_AL.csv",quote = F,row.names = F,sep=";")
write.table(my.KEGG.up@result,file = "KEGG_enrichment_AL.csv",quote = F,row.names = F,sep=";")
# my.Go.BP@result$Description[37]
# add filter, differen < -10, normal significant 
my.down.corrdinate<-intersect(which(x.new$difference< (-10)),which(abs(x.new$dist.to.feature)<=1000))
my.gene.down<-unique(x.new$Gene[my.down.corrdinate])
my.Go.CC.down<- enrichGO(gene=my.gene.down,OrgDb = org.Hs.eg.db,
                    ont="CC",keyType = "SYMBOL",pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

my.Go.BP.down<- enrichGO(gene=my.gene.down,OrgDb = org.Hs.eg.db,
                    ont="BP",keyType = "SYMBOL",pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

my.Go.MF.down<- enrichGO(gene=my.gene.down,OrgDb = org.Hs.eg.db,
                    ont="MF",keyType = "SYMBOL",pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
my.entrize.id.down<-AnnotationDbi::select(org.Hs.eg.db,keys = my.gene.down,
                                        keytype = "SYMBOL",
                                        columns = "ENTREZID")
my.KEGG.down<-enrichKEGG(gene=my.entrize.id.down$ENTREZID,organism = "hsa",
                       keyType = "kegg",pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
write.table(my.Go.CC.down@result,file = "GO_enrichment_CC_Normal.csv",quote = F,row.names = F,sep = ";")
write.table(my.Go.MF.down@result,file = "GO_enrichment_MF_Normal.csv",quote = F,row.names = F,sep=";")
write.table(my.Go.BP.down@result,file = "GO_enrichment_BP_Normal.csv",quote = F,row.names = F,sep=";")
write.table(my.KEGG.down@result,file = "KEGG_enrichment_Normal.csv",quote = F,row.names = F,sep=";")
# visualized in ggplot2
PlotPathway<-function(table1=NULL,table2=NULL,condition=c("event","normal"),cut.padj=0.05,PathwayType="KEGG",top=50,...){
  my.datafra1<-as.data.frame(table1@result)
  my.datafra2<-as.data.frame(table2@result)
  my.datafra1<-my.datafra1[my.datafra1$p.adjust<=cut.padj,]
  my.datafra1$condition<-condition[1]
  my.datafra2<-my.datafra2[my.datafra2$p.adjust<=cut.padj,]
  my.datafra2$condition<-condition[2]
  my.merged.datafra<-rbind.data.frame(my.datafra1,my.datafra2)
  my.merged.datafra<-my.merged.datafra[order(my.merged.datafra$p.adjust),]
  my.merged.datafra<-my.merged.datafra[1:top,]
  my.merged.datafra<-na.omit(my.merged.datafra)
  p<-ggplot(data=my.merged.datafra,aes(x=condition,y=Description))
  p<-p+geom_point(aes(size=my.merged.datafra$Count,color=my.merged.datafra$p.adjust))
  p<-p+scale_color_gradient(low = "black",high = "red",trans="reverse")
  p<-p+ylab("Pathway Discription")+labs(color="P-adj",size="counts",title = paste0(PathwayType))
  suppressWarnings(p)
}

PlotPathway(table1 = my.Go.BP.up,table2 = my.Go.BP.down,condition=c("Alzheimer","Normal"),top = 50,PathwayType = "GO_BP")
PlotPathway(table1 = my.Go.CC.up,table2 = my.Go.CC.down,condition=c("Alzheimer","Normal"),top = 50,PathwayType = "GO_CC")
PlotPathway(table1 = my.Go.MF.up,table2 = my.Go.MF.down,condition=c("Alzheimer","Normal"),top = 50,PathwayType = "GO_MF")
PlotPathway(table1 = my.KEGG.up,table2 = my.KEGG.down,condition=c("Alzheimer","Normal"),top = 50,PathwayType = "KEGG")







