# GO pathway 
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
### change DM_snps. csv gene name

gene.refseq.name<-unlist(sapply(strsplit(my.Diff.final$feature.name,"[.]"),'[[',1))
gene.Symbol<-AnnotationDbi::select(org.Hs.eg.db,keys = gene.refseq.name,
                                   keytype = "REFSEQ",
                                   columns = "SYMBOL")
x.new<-cbind.data.frame(my.Diff.final,feature.name=gene.refseq.name,Gene=gene.Symbol$SYMBOL)
write.csv(x.new,file = "DM_SNP_all.csv",quote = F,row.names = F)
x.new<-read.csv("DM_SNP_all.csv",header = T,stringsAsFactors = F)
# add filter qvalue<0.05
x.new<-x.new[x.new$qvalue<=0.05,]
# select feature.
# up regulate
# add filter, differen > 10
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
write.table(my.Go.CC@result,file = "GO_enrichment_CC_up.csv",quote = F,row.names = F,sep = ";")
write.table(my.Go.MF@result,file = "GO_enrichment_MF_up.csv",quote = F,row.names = F,sep=";")
write.table(my.Go.BP@result,file = "GO_enrichment_BP_up.csv",quote = F,row.names = F,sep=";")
# my.Go.BP@result$Description[37]
# down regulate
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
write.table(my.Go.CC@result,file = "GO_enrichment_CC_down.csv",quote = F,row.names = F,sep = ";")
write.table(my.Go.MF@result,file = "GO_enrichment_MF_down.csv",quote = F,row.names = F,sep=";")
write.table(my.Go.BP@result,file = "GO_enrichment_BP_down.csv",quote = F,row.names = F,sep=";")
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

PlotPathway(table1 = my.Go.BP.up,table2 = my.Go.BP.down,condition=c("Alzheimer","Normal"),top = 100,PathwayType = "GO_BP")










