library(Seurat)
library(dplyr)
library(ggplot2)
library(devtools)
library(knitr)
library(RColorBrewer)
library(lmerTest)

##Load object extracted from merged object
Monocytes <- readRDS('Monocytes.rds')


## Seurat analysis
# Normalize object
Monocytes <- NormalizeData(Monocytes)

# Scale object using all genes
all.genes<-rownames(Monocytes)
Monocytes <- ScaleData(Monocytes, features = all.genes)

## Find object variable genes 
Monocytes<-FindVariableFeatures(Monocytes, selection.method = "vst", nfeatures = 2000)

# Run PCA
Monocytes <- RunPCA(Monocytes,npcs = 80, ndims.print = 1:5)

# Use elbow plot to define number of dimensions
ElbowPlot(Monocytes, ndims = 80)

# run UMAP
Monocytes <- RunUMAP(object = Monocytes, dims = 1:80)

# plot basic UMAPs
p1<-DimPlot(object = Monocytes, reduction = "umap",pt.size = 0.1,label = TRUE, group.by = "orig.ident")

p2<-DimPlot(object = Monocytes, reduction = "umap",pt.size = 0.1,label = TRUE, group.by = "Patient") 

color<-brewer.pal(n = 11, name = "BrBG")
p3<-DimPlot(object = Monocytes, reduction = 'umap', label = FALSE, group.by ="Timepoint", pt.size = 0.1, 
            order = c('TA', 'TA2', 'TB', 'TC', 'TE', 'rec'),
            cols=rev(color))

p4<-DimPlot(object = Monocytes, reduction = "umap",pt.size = 0.1, label = F, group.by = "Pseudotime",
            order=c( "Healthy",'Recovered',"Complicated (incremental)",  "Complicated with hyperinflammatory syndrom",
                     "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic", 'Critical'),
            cols = c("#7031B9", "#FDC077", "#F99B1C", "#E65826",'#A04E9E',"#E84F8C",'#51BBFE','#A7A9AC'))
CombinePlots(plots = list(p1,p2,p3,p4), ncol=1)



## Clustering
# Find neighbours
Monocytes <- FindNeighbors(Monocytes,  dims = 1:80)

# Calculate clusters
Monocytes <- FindClusters(Monocytes, resolution = 0.1)

# Identify signature genes of each cluster
All_pre.markers <- FindAllMarkers(Monocytes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- All_pre.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
utop10<-unique(top10$gene)
color<-brewer.pal(n = 9, name = "RdBu")

# Plot top10 signature genes of each cluster
DotPlot(Monocytes, features=utop10,  dot.scale = 8) +
  scale_color_gradientn(colours  = rev(color)) + coord_flip() + scale_y_discrete(position = "right") +
  theme( axis.text.x = element_text(angle = 60, hjust = 0))


## Biomarkers of monocyte subtypes
# Cd14high vs. Cd16high
VlnPlot(object = Monocytes, features = c( "CD14", 'FCGR3A', 'ITGA4'), pt.size = 0)

# Classical monocytes
VlnPlot(object = Monocytes, features = c( "GPX1", "ATP5E", "GNB2L1", 'ATP5G2', 'ATP5L', 'ATP5I'), pt.size = 0)

# HLA-DRhigh CD83high
VlnPlot(object = Monocytes, features = c("IFI27",  "IFI44L", "PARP14", "CCNL1", 'NEAT1','RNF213'), pt.size = 0)

# HLA-DRhigh HBBhigh
VlnPlot(object = Monocytes, features = c('HBB',  'PF4'), pt.size = 0)

# CD163high
VlnPlot(object = Monocytes, features = c('ISG15', 'IFITM3', 'IFI6' , 'IFIT1', 'TNFSF10','MX1','IFIT3','LY6E'), pt.size = 0)

# S100high
VlnPlot(object = Monocytes, features = c('S100A12',  'CLU', 'CXCL8','MAFB','PLBD1','RGCC','RETN'), pt.size = 0)

# No-classical  monocytes
VlnPlot(object = Monocytes, features = c('FCGR3A',  'CDKN1C','MS4A7','LYPD2','NAP1L1','HES4'), pt.size = 0)


## Identification of monocyte subtypes 
# Add monocyte subtypes to individual clusters
Idents(Monocytes)<-Monocytes$seurat_clusters
new.cluster.ids <- c( 'Classical monocytes','Classical monocytes','Classical monocytes',
                      'Non-classical monocytes','CD163+ monocytes',
                      'S100A+ monocytes','Classical monocytes','S100A+ monocytes',
                      'Classical monocytes','CD163+ monocytes', 
                      'HLA-DR+ CD83+ monocytes', '?','HLA-DR+ HBB+ monocytes',
                      'CD163+ monocytes','Classical monocytes',
                      'S100A+ monocytes','Non-classical monocytes','Classical monocytes','double-positive monocytes', 'Classical monocytes')
names(new.cluster.ids) <- levels(Monocytes)
Monocytes <- RenameIdents(Monocytes, new.cluster.ids)
Monocytes$Celltype_monocytes<-Idents(Monocytes)

# Remove cluster that had no clear identification
Idents(Monocytes)<-Monocytes$Celltype_monocytes
Monocytes2<-subset(Monocytes, idents=c('Classical monocytes',
                                       'Non-classical monocytes',
                                       'CD163+ monocytes',
                                       'S100A+ monocytes',
                                       'HLA-DR+ CD83+ monocytes', 
                                       'HLA-DR+ HBB+ monocytes',
                                       'double-positive monocytes'))

# Plot monocytes subtypes as UMAP
color<-brewer.pal(n = 11, name = "BrBG")
DimPlot(Monocytes2, reduction = "umap", group.by = 'Celltype_monocytes',label=FALSE,
        cols=c( '#003c30', '#bf812d', '#dfc27d', '#c7eae5', '#80cdc1', '#35978f', '#543005'), pt.size = 0.1)


## Cell proportion
table<-read.csv2('~/Documents/IKMB/scRNAseq/ProjectC/Github_scripts/Prop_Celltypes_github.csv')

# Relative monocyte subtypes proportions (percentage)
table$rel_CM<-(table$CM/table$Total.cells)*100
table$rel_NCM<-(table$NCM/table$Total.cells)*100
table$rel_CD163<-(table$CD163/table$Total.cells)*100
table$rel_S100A<-(table$S100A/table$Total.cells)*100
table$rel_HLACD83<-(table$HLACD83/table$Total.cells)*100
table$rel_HLAHBB<-(table$HLAHBB/table$Total.cells)*100
table$rel_double<-(table$double/table$Total.cells)*100

# Plot monocyte subtypes proportions per pseudotime
p1 <- ggplot(table, aes(x=Pseudotime, y=rel_CM,  group=Pseudotime, color=Pseudotime)) 
p1 <- p1 +  geom_boxplot()
p1 <- p1 + geom_point(aes(color = Pseudotime, group=Pseudotime),na.rm = TRUE, size=8)
p1<- p1 + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                             axis.text.x = element_text(size=14, color = "black"), 
                             axis.title = element_text(size = 20))
p1<- p1 + scale_x_discrete(limits=c("Healthy","Complicated (incremental)",'Critical',"Complicated with hyperinflammatory syndrom",
                                    "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic",'Recovered'))
p1 <- p1 + scale_color_manual(values = c('#E84F8C',"#E65826","#A04E9E", "#7031B9","#A7A9AC", "#F99B1C","#51BBFE", "#FDC077"))
p1 <- p1 + theme_bw() +theme(legend.position = "none",
                             text = element_text(size=16),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())
p1 <- p1 + ggtitle("Classical monocytes")

p2 <- ggplot(table, aes(x=Pseudotime, y=rel_NCM,  group=Pseudotime, color=Pseudotime)) 
p2 <- p2 +  geom_boxplot()
p2 <- p2 + geom_point(aes(color = Pseudotime, group=Pseudotime),na.rm = TRUE, size=8)
p2 <- p2 + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                              axis.text.x = element_text(size=14, color = "black"), 
                              axis.title = element_text(size = 20))
p2<- p2 + scale_x_discrete(limits=c("Healthy","Complicated (incremental)",'Critical',"Complicated with hyperinflammatory syndrom",
                                    "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic",'Recovered'))
p2 <- p2 + scale_color_manual(values = c('#E84F8C',"#E65826","#A04E9E", "#7031B9","#A7A9AC", "#F99B1C","#51BBFE", "#FDC077"))
p2 <- p2 + theme_bw() +theme(legend.position = "none",
                             text = element_text(size=16),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())
p2 <- p2 + ggtitle("Non-classical monocytes")

p3 <- ggplot(table, aes(x=Pseudotime, y=rel_CD163,  group=Pseudotime, color=Pseudotime)) 
p3 <- p3 +  geom_boxplot()
p3 <- p3 + geom_point(aes(color = Pseudotime, group=Pseudotime),na.rm = TRUE, size=8)
p3 <- p3 + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                              axis.text.x = element_text(size=14, color = "black"), 
                              axis.title = element_text(size = 20))
p3<- p3 + scale_x_discrete(limits=c("Healthy","Complicated (incremental)",'Critical',"Complicated with hyperinflammatory syndrom",
                                    "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic",'Recovered'))
p3 <- p3 + scale_color_manual(values = c('#E84F8C',"#E65826","#A04E9E", "#7031B9","#A7A9AC", "#F99B1C","#51BBFE", "#FDC077"))
p3 <- p3 + theme_bw() +theme(legend.position = "none",
                             text = element_text(size=16),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())
p3 <- p3 + ggtitle("CD163+ monocytes")

p4 <- ggplot(table, aes(x=Pseudotime, y=rel_S100A,  group=Pseudotime, color=Pseudotime)) 
p4 <- p4 +  geom_boxplot()
p4 <- p4 + geom_point(aes(color = Pseudotime, group=Pseudotime),na.rm = TRUE, size=8)
p4 <- p4 + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                              axis.text.x = element_text(size=14, color = "black"), 
                              axis.title = element_text(size = 20))
p4<- p4 + scale_x_discrete(limits=c("Healthy","Complicated (incremental)",'Critical',"Complicated with hyperinflammatory syndrom",
                                    "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic",'Recovered'))
p4 <- p4 + scale_color_manual(values = c('#E84F8C',"#E65826","#A04E9E", "#7031B9","#A7A9AC", "#F99B1C","#51BBFE", "#FDC077"))
p4 <- p4 + theme_bw() +theme(legend.position = "none",
                             text = element_text(size=16),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())
p4 <- p4 + ggtitle("S100A+  monocytes")

p5 <- ggplot(table, aes(x=Pseudotime, y=rel_HLACD83,  group=Pseudotime, color=Pseudotime)) 
p5 <- p5 +  geom_boxplot()
p5 <- p5 + geom_point(aes(color = Pseudotime, group=Pseudotime),na.rm = TRUE, size=8)
p5 <- p5 + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                              axis.text.x = element_text(size=14, color = "black"), 
                              axis.title = element_text(size = 20))
p5<- p5 + scale_x_discrete(limits=c("Healthy","Complicated (incremental)",'Critical',"Complicated with hyperinflammatory syndrom",
                                    "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic",'Recovered'))
p5 <- p5 + scale_color_manual(values = c('#E84F8C',"#E65826","#A04E9E", "#7031B9","#A7A9AC", "#F99B1C","#51BBFE", "#FDC077"))
p5 <- p5 + theme_bw() +theme(legend.position = "none",
                             text = element_text(size=16),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())
p5 <- p5 + ggtitle("HLA-DR+ CD83+ monocytes")

p6 <- ggplot(table, aes(x=Pseudotime, y=rel_CM,  group=Pseudotime, color=Pseudotime)) 
p6 <- p6 +  geom_boxplot()
p6 <- p6 + geom_point(aes(color = Pseudotime, group=Pseudotime),na.rm = TRUE, size=8)
p6 <- p6 + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                              axis.text.x = element_text(size=14, color = "black"), 
                              axis.title = element_text(size = 20))
p6<- p6 + scale_x_discrete(limits=c("Healthy","Complicated (incremental)",'Critical',"Complicated with hyperinflammatory syndrom",
                                    "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic",'Recovered'))
p6 <- p6 + scale_color_manual(values = c('#E84F8C',"#E65826","#A04E9E", "#7031B9","#A7A9AC", "#F99B1C","#51BBFE", "#FDC077"))
p6 <- p6 + theme_bw() +theme(legend.position = "none",
                             text = element_text(size=16),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())
p6 <- p6 + ggtitle("HLA-DR+ HBB+ monocytes")

p7 <- ggplot(table, aes(x=Pseudotime, y=rel_CM,  group=Pseudotime, color=Pseudotime)) 
p7 <- p7 +  geom_boxplot()
p7 <- p7 + geom_point(aes(color = Pseudotime, group=Pseudotime),na.rm = TRUE, size=8)
p7 <- p7 + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                              axis.text.x = element_text(size=14, color = "black"), 
                              axis.title = element_text(size = 20))
p7<- p7 + scale_x_discrete(limits=c("Healthy","Complicated (incremental)",'Critical',"Complicated with hyperinflammatory syndrom",
                                    "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic",'Recovered'))
p7 <- p7 + scale_color_manual(values = c('#E84F8C',"#E65826","#A04E9E", "#7031B9","#A7A9AC", "#F99B1C","#51BBFE", "#FDC077"))
p7 <- p7 + theme_bw() +theme(legend.position = "none",
                             text = element_text(size=16),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())
p7<- p7 + ggtitle("CD14+ CD16+ monocytes")

CombinePlots(plots = list(p1,p2,p3,p4,p5, p6, p7), ncol=1)


## Different cell proportions through pseudotimes
# select patients with covid19
table2<-subset(table, table$Pseudotime %in% c("Complicated (incremental)",'Critical',"Complicated with hyperinflammatory syndrom",
                                              "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic",'Recovered'))

# Longitudinal liner mixed model
fit <- lmer(rel_CM ~ Pseudotime + (1|Patient), data=table2, REML = FALSE)
reduced.fit<- lmer(rel_CM ~ 1 + (1|Patient),  data=table2, REML = FALSE)
anova(reduced.fit,fit )

# select control patients
table_healthy<-subset(table, table$Pseudotime %in% c('Healthy'))

#Mann-whitney test comparing healthy with covid19 patients
wilcox.test(table2$rel_double, table_healthy$rel_double)


## Signature genes per pseudotime for monocytes subtypes of interest

#Subset monocytes subtypes of interest
Idents(Monocytes)<-Monocytes$Celltype_monocytes
CM<-subset(Monocytes, idents='Classical monocytes')
NCM<-subset(Monocytes, idents='Non-classical monocytes')

# sort object by pseudotime
Idents(CM) <- CM$Pseudotime
my_levels <- c("Healthy","Complicated (incremental)",'Critical', "Complicated with hyperinflammatory syndrom",
               "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic", 'Recovered')
Idents(CM) <- factor(Idents(CM), levels= my_levels)
Idents(NCM) <- NCM$Pseudotime
my_levels <- c("Healthy","Complicated (incremental)",'Critical', "Complicated with hyperinflammatory syndrom",
               "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic", 'Recovered')
Idents(NCM) <- factor(Idents(NCM), levels= my_levels)

# Calculate signature genes per pseudotime for Classical monocytes
All_pre.markers <- FindAllMarkers(CM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- All_pre.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
utop10<-unique(top10$gene)
color<-brewer.pal(n = 9, name = "RdBu")

# Plot top10 signature genes per pseudotime for Classical monocytes
DotPlot(CM, features=utop10,  dot.scale = 8) +
  scale_color_gradientn(colours  = rev(color)) + coord_flip() + scale_y_discrete(position = "right") +
  theme( axis.text.x = element_text(angle = 60, hjust = 0))

# Calculate signature genes per pseudotime for Non-classical monocytes
All_pre.markers <- FindAllMarkers(CM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- All_pre.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
utop10<-unique(top10$gene)
color<-brewer.pal(n = 9, name = "RdBu")

# Plot top10 signature genes per pseudotime for Non-classical monocytes
DotPlot(CM, features=utop10,  dot.scale = 8) +
  scale_color_gradientn(colours  = rev(color)) + coord_flip() + scale_y_discrete(position = "right") +
  theme( axis.text.x = element_text(angle = 60, hjust = 0))


#Plot HLA genes of interest for Classical monocytes
Idents(CM)<-CM$Pseudotime
my_levels <- c("Healthy","Complicated (incremental)",'Critical', "Complicated with hyperinflammatory syndrom",
               "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic", 'Recovered')
Idents(CM) <- factor(Idents(CM), levels= my_levels)

VlnPlot(CM, features =c('HLA-DRA','HLA-DRB5', 'HLA-DQA1', 'HLA-DQA2',
                        'HLA-DPA1', 'HLA-DPB1'),  pt.size = 0,ncol = 3,
        cols = c("#A7A9AC", '#E84F8C', "#7031B9", "#A04E9E",'#E65826',"#F99B1C",'#FDC077','#51BBFE'))

#Plot HLA genes of interest for Non-classical monocytes

Idents(NCM)<-NCM$Pseudotime
my_levels <- c("Healthy","Complicated (incremental)",'Critical', "Complicated with hyperinflammatory syndrom",
               "Complicated (recovering)","Mild (recovering)","recovery/asymptomatic", 'Recovered')
Idents(NCM) <- factor(Idents(NCM), levels= my_levels)

VlnPlot(NCM, features =c('HLA-DRA','HLA-DRB5', 'HLA-DQA1', 'HLA-DQA2',
                         'HLA-DPA1', 'HLA-DPB1'),  pt.size = 0,ncol = 3,
        cols = c("#A7A9AC", '#E84F8C', "#7031B9", "#A04E9E",'#E65826',"#F99B1C",'#FDC077','#51BBFE'))

