wkdir <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path))

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",)
BiocManager::install("DESeq2","GO.db")
BiocManager::install("GO.db")
BiocManager::install("WGCNA")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
library(DESeq2)
library(WGCNA)

########Figure 4#############
FC_count<-read.csv("../Data/FC_count.csv")
coldata<-read.csv("../Data/coldata.csv", head=T,row.names=1,stringsAsFactors = F)
Blue<-read.csv("../Data/Blue.csv", head=T,stringsAsFactors = F)
Green<-read.csv("../Data/Green.csv", head=T,stringsAsFactors = F)
Purple<-read.csv("../Data/Purple.csv", head=T,stringsAsFactors = F)
Metabolite_serum<-read.csv("../Data/Metabolite_serum.csv")

#preparing count file from Featurecounts
rownames(FC_count)<-FC_count$gene_id
FC_count$gene_id<-NULL
FC_count<-as.matrix(FC_count)

#preparing coldata
rownames(coldata)<-coldata$ID
coldata$ID<-NULL
coldata$group<-paste0(coldata$Treatment,sep="_D",coldata$Day)
coldata$group<-factor(coldata$group,levels=c("Control_D0","Infected_D2","Infected_D4","Infected_D7","Infected_D14"))
coldata$Treatment<-factor(coldata$Treatment,levels=c("Control","Infected"))
coldata$Day<-factor(coldata$Day,levels=c("0",'2','4','7','14'))

#DESEQ2 analysis
dds<-DESeqDataSetFromMatrix(countData=FC_count,colData = coldata, design=~group)
dds<-DESeq(dds)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

##########Figure 4b############
plotPCA(vsd,intgroup="group")

#obtain DEGs
res.D2<-results(dds,contrast=c('group','Infected_D2','Control_D0'),alpha=0.05)
res.D2<-res.D2[order(res.D2$padj),]
sig.res.D2<-res.D2[which((res.D2$padj<0.05) & (abs(res.D2$log2FoldChange)>1)),]

res.D4<-results(dds,contrast=c('group','Infected_D4','Control_D0'),alpha=0.05)
res.D4<-res.D4[order(res.D4$padj),]
sig.res.D4<-res.D4[which((res.D4$padj<0.05) & (abs(res.D4$log2FoldChange)>1)),]

res.D7<-results(dds,contrast=c('group','Infected_D7','Control_D0'),alpha=0.05)
res.D7<-res.D7[order(res.D7$padj),]
sig.res.D7<-res.D7[which((res.D7$padj<0.05) & (abs(res.D7$log2FoldChange)>1)),]

res.D14<-results(dds,contrast=c('group','Infected_D14','Control_D0'),alpha=0.05)
res.D14<-res.D14[order(res.D14$padj),]
sig.res.D14<-res.D14[which((res.D14$padj<0.05) & (abs(res.D14$log2FoldChange)>1)),]

##########Figure 4C#############

library('RColorBrewer')
library('pheatmap')

DEGlist<-union(rownames(sig.res.D2),row.names(sig.res.D4))
DEGlist<-union(DEGlist,row.names(sig.res.D7))
DEGlist<-union(DEGlist,row.names(sig.res.D14))

vsdtable<-as.data.frame(assay(vsd))
vsdtable<-vsdtable[,c(1:5,11:24,6:10)]
vsdtable_DEGlist<-subset(vsdtable,rownames(vsdtable) %in% DEGlist)
vsdmat<-as.matrix(vsdtable_DEGlist)
write.csv(vsdmat, "../Data/vsdmat.csv")

coldata_hm<-coldata[c(1:5,11:24,6:10),]
coldata_hm<-data.frame(coldata_hm[,c(3),drop=FALSE])

my_color_annotation<-list(group = c(Control_D0 = "#FFFF08", 
                                    Infected_D2 = "#FFDA08",
                                    Infected_D4 ="#FFB508",
                                    Infected_D7 ="#FF8C08",Infected_D14 ="#FF5508"),
                          cluster=c('5'="#99E600",'1'="#68DC00",'2'="#00D43D",
                                    '7'="#00978C",'3'="#095E9F",'6'="#132AA9",'4'="#6D08A5"))
hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

HM_vsdmat<-pheatmap(vsdmat, scale = 'row', border_color = 'NA',
                    color = hm_color,show_rownames = FALSE,
                    cluster_rows = TRUE, cluster_cols = FALSE,
                    annotation_col=coldata_hm,
                    annotation_colors =my_color_annotation,
                    clustering_distance_rows = 'correlation',
                    angle_col =45)

clusterlist<-as.data.frame(vsdmat[HM_vsdmat$tree_row[['order']],]) #obtain gene clusters #each gene cluster is exported for functional enrichment analysis
clus<-as.data.frame(cutree(HM_vsdmat$tree_row, k=7))
colnames(clus)<-c('cluster')
clus<-clus[rownames(clusterlist),,drop=FALSE]
clusterlist$cluster<-clus$cluster
roldata_hm<-clusterlist[,c(25),drop=FALSE]
roldata_hm$cluster<-factor(roldata_hm$cluster,levels = c("5",'1','2','7','3','6','4'))

HM_vsdmat<-pheatmap(vsdmat, scale = 'row', border_color = 'NA',
                    color = hm_color,show_rownames = FALSE,
                    cluster_rows = TRUE, cluster_cols = FALSE,
                    annotation_col=coldata_hm,
                    annotation_row = roldata_hm,
                    annotation_colors = my_color_annotation,
                    clustering_distance_rows = 'correlation',
                    angle_col =45)

#########Figure 6a################
#serum metabolite matrix
Metabolite_serum<-read.csv("../Data/Metabolite_serum.csv",head=T,stringsAsFactors = F)

rownames(Metabolite_serum)<-Metabolite_serum$ID
Metabolite_serum<-Metabolite_serum[,-c(1:4)]

hm_color_cor<- colorRampPalette(c("blue", "white", "red"))(100)

HM_Metabolite_serum<-pheatmap(Metabolite_serum, scale = 'column', border_color = 'NA',
                         color = hm_color_cor,show_rownames = TRUE,cellwidth = 12,
                         cellheight = 12,
                         cluster_rows = FALSE, cluster_cols = TRUE,
                         angle_col =45,cutree_cols = 2)

#########Figure 6b################
#metabolite group data #sum of abundance 

Metabolite_group<-as.data.frame(Metabolite_serum[,HM_Metabolite_serum$tree_col[['order']]])
Metabolite_group<-as.data.frame(t(Metabolite_group))
  
group<-as.data.frame(cutree(HM_Metabolite_serum$tree_col, k=2))
colnames(group)<-c('group')
#group<-group[colnames(Metabolite_group),,drop=FALSE]

Metabolite_group$group<-group$group

Metabolite_group1<-subset(Metabolite_group,Metabolite_group$group=="1")
Metabolite_group1$group<-NULL
Metabolite_group1<-as.data.frame(t(Metabolite_group1))
Metabolite_group1$Group1<-rowSums(Metabolite_group1)
Metabolite_group1<-Metabolite_group1[,c(19),drop=F]

Metabolite_group2<-subset(Metabolite_group,Metabolite_group$group=="2")
Metabolite_group2$group<-NULL
Metabolite_group2<-as.data.frame(t(Metabolite_group2))
Metabolite_group2$Group2<-rowSums(Metabolite_group2)
Metabolite_group2<-Metabolite_group2[,c(14),drop=F]

Metabolite_abundance<-as.data.frame(cbind(Metabolite_group1,Metabolite_group2))

###WGCNA######
library(WGCNA)

exprsdat<vsdmat[,c(1:19)]
exprsdat<-exprsdat[,-c(10)] #ensure the samples are matched with those from serum metabolomics
exprsdat<-as.data.frame(t(exprsdat))
##check missing values
###if result is TRUE, then all genes have passed the test
gsg = goodSamplesGenes(exprsdat, verbose = 3);
gsg$allOK

#omit unqualified gene
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(exprsdat)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(exprsdat)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exprsdat = exprsdat[gsg$goodSamples, gsg$goodGenes]
}


TranscriptomeSamples = rownames(exprsdat)
traitRows = match(TranscriptomeSamples, rownames(Metabolite_abundance))
datTraits = Metabolite_abundance[traitRows, c(1:2)]

clusterlist2<-clusterlist[colnames(exprsdat),]
moduleColors = clusterlist2$cluster

# Calculate eigengenes of different clusters
MEs0 = moduleEigengenes(exprsdat, moduleColors)$eigengenes
MEs = MEs0
MEs = MEs[,rev(c(5,1,2,7,3,6,4))]
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Define numbers of genes and samples
nGenes = ncol(exprsdat);
nSamples = nrow(exprsdat);

sizeGrWindow(6,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot

Cor_HM<-labeledHeatmap(Matrix = moduleTraitCor,
                       xLabels = names(datTraits),
                       yLabels = names(MEs),
                       ySymbols = names(MEs),
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = 0.5,
                       zlim = c(-1,1),
                       bg.lab.y = rev(c("#99E600","#68DC00","#00D43D","#00978C","#095E9F","#132AA9","#6D08A5")),
                       main = paste("Module-trait relationships"))


##########Figure 6C#############
library(ggplot2)
library(RColorBrewer)

colnames(Green)<-c("cat",'logp','ratio','zscore','molecules')
Green$percent<-Green$ratio*100
Green<-Green[,-c(4:5)]
Green$clus<-rep("Green",10)

colnames(Blue)<-c("cat",'logp','ratio','zscore','molecules')
Blue$percent<-Blue$ratio*100
Blue<-Blue[,-c(4:5)]
Blue$clus<-rep("Blue",10)

colnames(Purple)<-c("cat",'logp','ratio','zscore','molecules')
Purple$percent<-Purple$ratio*100
Purple<-Purple[,-c(4:5)]
Purple$clus<-rep("Purple",10)

IPA<-rbind(Green,Blue,Purple)
IPA$clus<-factor(IPA$clus,levels=c("Green",'Blue','Purple'))

IPA_plot<-ggplot(data=IPA, aes(x=percent,y=reorder(cat,percent),color=logp,size=percent))+
  geom_point()+
  xlim(c(0,30))+
  scale_color_gradientn(colours=colorRampPalette(c("#4951a3",'red'))(20),
                        breaks=c(2,3,4,5),labels=c(2,3,4,5),limits=c(1,6))+
  facet_grid(rows=vars(clus),scales = "free_y")+
  scale_size(range = c(4,10))+
  theme_gray()+
  labs(x="Percentage enrichment(%)",y="IPA canonical pathway",color="-log (p value)",size="Percentage enrichment(%)")
IPA_plot
