library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ImpulseDE2)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)


pseudotime_colors <- c("#A7A9AC","#E84F8C", "#7031B9", "#A04E9E", "#E65826", "#F99B1C", "#FDC077", "#51BBFE")

#Code for generating PCA plot using transcriptome data---------

##Read sample info file
covid_info_file <- read.csv("../info_files/Sequenced_sample_info_020720.csv")
covid_info_file <- subset(covid_info_file, covid_info_file$RNAseqID != '')

##Read raw count file
count_data <- read.csv("count_files/COVID_control_raw_counts_070720.txt", sep = '\t')
count_data <- count_data[, as.character(covid_info_file$RNAseqID)]

##Generate sample col data
col_data <- covid_info_file
rownames(col_data) <- col_data$RNAseqID
col_data$Gender <- as.factor(col_data$Gender)
col_data$Disease.traj <- as.factor(col_data$Disease.traj)
col_data$Pseudotime <- as.factor(as.character(col_data$Pseudotime))

##Create DESeq object
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Pseudotime)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

##Perform rlog transformation
rlog_counts <- rlog(dds_counts, blind = TRUE)
rlog.norm.counts <- assay(rlog_counts)

##Perform principle component analysis
pc <-prcomp(t(rlog.norm.counts))
pc_df <- as.data.frame(pc$x)
pc_df$Pseudotime <- col_data$Pseudotime
pc_df$Patient <- col_data$PatientID

summary(pc)

##Code to generate figure 3B
pdf("Control_COVID_PCA_plot.pdf", width = 10, height = 8)
P <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC2, color=Pseudotime, label=Patient)) + 
  geom_point(size=5) + geom_text(hjust=0, nudge_x = 2, nudge_y = 1, cex=3) + xlab("PC1 (20.35% variance)") + 
  ylab("PC2 (14.81% variance)")
P <- P + scale_color_manual(values = pseudotime_colors)
P <- P + theme_bw() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22))
P
dev.off()

#Code for pairwise differential expression analysis------

#Read reference gene info file
mart_export <- read.table("../..//reference/human/gencode_v25_gene_info.txt", header = T, sep="\t")
unique_mart <- subset(mart_export, duplicated(mart_export$Gene_id) == FALSE)
rownames(unique_mart) <- unique_mart$Gene_id

pseudotime <- c("1", "2", "3", "4", "5", "6", "7")
folder <- "output/DESeq2_results/Control_vs_COVID"
if (!dir.exists(folder)){
  dir.create(folder)
}

#Compare each pseudotime with control samples
for (i in 1:length(pseudotime)) {
  output_folder <- file.path(folder, paste("Control_vs_COVID", pseudotime[i], sep = ''))
  if (!dir.exists(output_folder)){
    dir.create(output_folder)
  }
  
  #Filter info data for desired pseudotime
  col_data <- subset(covid_info_file, covid_info_file$Pseudotime %in% c('0', pseudotime[i]))
  rownames(col_data) <- col_data$RNAseqID
  col_data$Pseudotime <- as.factor(as.character(col_data$Pseudotime))
  
  #Filter count data for desired pseudotime
  selected_count_data <- count_data[, as.character(col_data$RNAseqID)]
  
  #Create deseq object
  dds_counts <- DESeqDataSetFromMatrix(countData = selected_count_data, colData = col_data, design = ~ Pseudotime)
  dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
  dds_counts <- estimateSizeFactors(dds_counts)
  
  #Run DESeq
  dds <- DESeq(dds_counts, betaPrior = FALSE)
  res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
  res_sorted <- res[order(res$padj), ]
  
  #Add gene symbols
  res_sorted$gene <- unique_mart[rownames(res_sorted), ]$Gene_name
  write.table(res_sorted, file = file.path(output_folder, paste("DESeq2result_genenames_Control_vs_COVID_", pseudotime[i], ".txt", sep = '')), 
              sep="\t", quote=FALSE)
  top_gene <- rownames(res_sorted)[1]
}




#Code to generate volcano plot as in figure 3C------------

plot_data <- data.frame()

#Read each deseq out file and filter for significant differentially expressed genes
for (i in 1:length(pseudotime)) {
  input_folder <- file.path(folder, paste("Control_vs_COVID", pseudotime[i], sep = ''))
  result <- read.csv(file.path(input_folder, paste("DESeq2result_genenames_Control_vs_COVID_", pseudotime[i], ".txt", sep = '')), 
                     sep = '\t')
  result <- subset(result, result$padj < 0.05 & result$baseMean > 100 & abs(result$log2FoldChange) > 0.5)
  result$pseudotime <- pseudotime[i]
  result <- rownames_to_column(result, "Gene_ID")
  plot_data <- rbind(plot_data, result)
  
}

#Add direction of regulation to plot_data
plot_data$direction <- ifelse(plot_data$log2FoldChange < 0, "down", "up")

#Calculate the number of comparisons in which a gene is significant
gene_frequency <- table(plot_data$Gene_ID)
plot_data$frequency <- as.vector(gene_frequency[as.character(plot_data$Gene_ID)])

#Select interesting genes to label
gene_to_label <- c("MAOA", "HTRA3", "ERG", "TMEM119", "INHBB", "PTX3", "FCGR1CP", "PDZK1IP1", "FHL2", "NRIR", "SLC4A1", "C1QB", "TESC", "MERTK", 
                   "IDO1", "P2RY14", "F13A1", "STAT1", "GATA2", "NFKBIA", "TREML1", "NCR1", "NLRP6", "NLRP3", "JAK3", "TXK", "IL5RA", "SATB1-AS1", 
                   "SHARPIN", "TNFRSF10C", "CXCR1", "ATF6", "IFI16", "AHR", "IL13RA1", "C5AR1", "CSF2RA", "FYN", "IGHG1", "JCHAIN", "TOP2A", 
                   "IFI27", "BIRC5", "IGHA1", "IL1R2", "DEFA4", "DEFA3", "LTF", "RORC", "IFIT1B", "AIM2", "IL24", "CD38", "TNFRSF13C", "ERBB2", 
                   "TCF7", "CIITA", "TGFBI", "SIGIRR", "IGHA2", "CDC6", "CHI3L1", "HLA-DPA1", "HLA-DMB", "CD4", "IL6ST", "CD274", "C3AR1", "ADAM17", 
                   "C1QA", "CARD11", "CD14", "BCL3", "TICAM1")
gene_to_label2 <- c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DPB1", "HLA-DPB", "HLA-DQA1", "HLA-DQB1", "HLA-DR1", 
                    "HLA-DRB1", "HLA-DRB5")

#Too few DEGs at pseudotime 7, hence removed
plot_data1 <- subset(plot_data, plot_data$pseudotime != "7")

#Code to plot 
cairo_ps(file = "deg_pseudotime_lfc.ps", width = 9, height = 7)
p <- ggplot(data = plot_data1, aes(x=pseudotime, y=log2FoldChange, fill = direction))
p <- p + geom_point(position=position_jitter(width=0.1), aes(size=-1*log10(padj), alpha=frequency), pch=21, color="black", stroke=0.1)
p <- p + geom_text_repel(data = subset(plot_data1, plot_data1$gene %in% gene_to_label & abs(plot_data1$log2FoldChange) > 2), aes(x=pseudotime, y=log2FoldChange, label=gene), 
                         color="black", size=3.5, nudge_x = -0.5)
p <- p + geom_text_repel(data = subset(plot_data1, plot_data1$gene %in% gene_to_label2& abs(plot_data1$log2FoldChange) > 1), aes(x=pseudotime, y=log2FoldChange, label=gene), 
                         color="black", size=3.5, nudge_x = -0.5)
p <- p + scale_color_manual(values=c('#004C99', '#990000')) + scale_fill_manual(values=c('#004C99', '#990000'))
p <- p + xlab("Pseudotime") + ylab("Log fold change") + scale_x_discrete(labels = c("0 vs 1", "0 vs 2", "0 vs 3", "0 vs 4", "0 vs 5", "0 vs 6"))
p <- p + theme_bw() + theme(axis.text=element_text(size=14, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=12), 
                            plot.title = element_text(size=20, face="bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
p
dev.off()




#Code for longitudinal differential expression analysis---------

#Include only longitudinal samples
col_data <- subset(covid_info_file, !(covid_info_file$Pseudotime %in% c(0, 2, 7)))
col_data$Pseudotime <- as.factor(as.character(col_data$Pseudotime))

#Filter count data accordingly
covid_count_data <- count_data[, as.character(col_data$RNAseqID)]
colnames(covid_count_data) <- col_data$Name

#Create deseq object
dds_counts <- DESeqDataSetFromMatrix(countData = covid_count_data, colData = col_data, design = ~ PatientID + Pseudotime)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]


#Create impulseDE2 object and run
matCountData <- counts(dds_counts)
dfAnnotation <- data.frame(Sample=col_data$Name, Condition=rep('case', nrow(col_data)), 
                           Pseudotime=col_data$Pseudotime, Batch=col_data$PatientID, 
                           Time=as.numeric(as.character(col_data$Pseudotime)))
dfAnnotation$Time <- ifelse(dfAnnotation$Pseudotime == '1', 1, as.numeric(as.character(dfAnnotation$Pseudotime))-1)
dfAnnotation <- dfAnnotation[order(dfAnnotation$Time),]
matCountData <- matCountData[, as.character(dfAnnotation$Sample)]

objectImpulseDE2 <- runImpulseDE2(
  matCountData           = matCountData, 
  dfAnnotation           = dfAnnotation,
  boolCaseCtrl           = FALSE,
  vecConfounders         = c("Batch"),
  boolIdentifyTransients = TRUE,
  scaNProc               = 2 )

res_sorted <- objectImpulseDE2$dfImpulseDE2Results[order(objectImpulseDE2$dfImpulseDE2Results$padj),]
res_sorted$gene <- unique_mart[rownames(res_sorted), ]$Gene_name
write.table(res_sorted2, file = "output/Longitudinal_analysis/COVID_ImpulseDE2_result.txt", sep="\t", quote=FALSE)


#Code for comparison between pseudotime 2 and 3-------

#Filter info data for pseudotimes 2 and 3
col_data <- subset(covid_info_file, covid_info_file$Pseudotime %in% c(2,3))
col_data$Pseudotime <- as.factor(as.character(col_data$Pseudotime))

#Filter count data accordingly
covid_count_data <- count_data[, as.character(col_data$RNAseqID)]
colnames(covid_count_data) <- col_data$Name

#Create deseq object
dds_counts <- DESeqDataSetFromMatrix(countData = covid_count_data, colData = col_data, design = ~ Pseudotime)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)


#Run DESeq
dds <- DESeq(dds_counts, betaPrior = FALSE)
res <- results(dds, independentFiltering = TRUE, alpha = 0.05)

res_sorted$gene <- unique_mart[rownames(res_sorted), ]$Gene_name
write.table(res_sorted, file = "output//DESeq2_results/Pseudotime_2_VS_3/DESeq2result_genenames.txt", 
            sep="\t", quote=FALSE)

#Plot heatmap as in figure S2B
top_up <- subset(res_sorted, res_sorted$padj < 0.05 & res_sorted$log2FoldChange > 0)

top_down <- subset(res_sorted, res_sorted$padj < 0.05 & res_sorted$log2FoldChange < 0)

degs <- c(rownames(top_up), rownames(top_down))
matDataHeat <- normalized_counts[as.character(degs), ]
rownames(matDataHeat) <- c(as.character(top_up$gene), as.character(top_down$gene))

matDataHeat <- t(scale(t(matDataHeat)))

col_ha = HeatmapAnnotation(df = data.frame(Pseudotime = covid_info_file$Pseudotime),
                           col = list(Pseudotime = c('2' = '#7031B9', '3' = '#A04E9E')), 
                           gp = gpar(col = "black"))
pdf("output/DESeq2_results/Pseudotime_2_vs_3/DEGs_heatmap.pdf", width = 7, height = 9)
Heatmap(matDataHeat, col=rev(brewer.pal(11, "RdYlBu")), top_annotation = col_ha)
dev.off()

#Code to generate heatmap as in figure 3D-------
covid_info_file$Sample_name <- paste(covid_info_file$Name, covid_info_file$Pseudotime, sep = '_')

##Read normalized read counts
normalized_expression_counts <- read.csv("count_files/COVID_control_normalized_counts_070720.txt", sep = '\t')
normalized_expression_counts <- normalized_expression_counts[, as.character(covid_info_file$RNAseqID)]
colnames(normalized_expression_counts) <- covid_info_file$Sample_name

#Get list of longitudinal DEGs
impulse_degs <- read.csv("output/Longitudinal_analysis/COVID_ImpulseDE2_result.txt", 
                         sep = '\t')
impulse_degs <- subset(impulse_degs, impulse_degs$padj < 0.05)
if (nrow(impulse_degs) > 500) {
  impulse_degs <- impulse_degs[1:500,]
}
impulse_degs <- rownames(impulse_degs)
severe_degs <- read.csv("output/DESeq2_results/Pseudotime_2_vs_3/DESeq2result_genenames.txt", sep = '\t')
severe_degs <- subset(severe_degs, severe_degs$padj < 0.05)
severe_degs <- rownames(severe_degs)

sig_genes <- union(impulse_degs, severe_degs)

##Get expression data for DEGs
sig_data <- normalized_expression_counts[sig_genes, ]

genenames <- as.character(unique_mart[sig_genes, ]$Gene_name)

##Average expression values for each pseudotime
vecTP <- str_split_fixed(colnames(sig_data), '_', 4)[, 4]
vecUniqueTP <- unique(vecTP)

matDataHeat <- do.call(cbind, lapply(vecUniqueTP, function(tp) {
  vecidxCols <- which(vecTP %in% tp)
  if (length(vecidxCols) > 1) {
    return(rowMeans(sig_data[, vecidxCols], na.rm = TRUE))
  } else {
    return(sig_data[, vecidxCols])
  }
}))
colnames(matDataHeat) <- vecUniqueTP
rownames(matDataHeat) <- genenames

#Generate heatmap
matDataHeat <- t(scale(t(matDataHeat)))
matDataHeat <- matDataHeat[, c("0", "1", "3", "4", "5", "6", "2")]

col_ha = HeatmapAnnotation(df = data.frame(Pseudotime = colnames(matDataHeat)),
                           col = list(Pseudotime = c('0' = 'grey', '1' = '#E84F8C', '2' = '#7031B9', '3' = '#A04E9E', 
                                                     '4' = '#E65826', '5' = '#F99B1C', '6' = '#FDC077')), 
                           gp = gpar(col = "black"))

pdf("output/impulse+severe_degs_heatmap.pdf", width=6, height=12)
Heatmap(matDataHeat, col=rev(brewer.pal(11, "RdYlBu")), bottom_annotation = col_ha, cluster_columns = FALSE, column_split = c("A", rep("B", 5), "C"), 
        clustering_distance_rows = function(x, y) 1 - cor(x, y), row_km = 7, row_gap = unit(1, "mm"), row_names_gp = gpar(fontsize = 2), border = TRUE)
dev.off()

#Code to generate bar plot with number of DEGs as in figure S2A----
deg_data <- data.frame(timepoint=as.character(), upgenes=as.integer(), downgenes=as.integer())

##Get number of DEGs at each pseudotime compared to controls
for (i in 1:length(pseudotime)) {
  input_folder <-  paste("output/DESeq2_results/Control_vs_COVID/Control_vs_COVID", pseudotime[i], sep = '')
  result <- read.csv(file.path(input_folder, paste("DESeq2result_genenames_Control_vs_COVID_", pseudotime[i], ".txt", sep = '')), 
                     sep = '\t')
  result <- subset(result, result$padj < 0.05 & result$baseMean > 100 & abs(result$log2FoldChange) > 0.5)
  uplen <- nrow(subset(result, result$log2FoldChange > 0))
  downlen <- nrow(subset(result, result$log2FoldChange < 0))
  deg_data <- rbind(deg_data, data.frame(pseudotime[i], uplen, downlen))
}
##Reformat data
deg_data <- melt(deg_data)

##Get number of longitudinal DEGs
impulse_data <- read.csv("output/Longitudinal_analysis/COVID_ImpulseDE2_result.txt", 
                         header = TRUE, sep = '\t')
impulse_data <- subset(impulse_data, impulse_data$padj < 0.05)

deg_data <- rbind(deg_data, data.frame(pseudotime.i.="Longitudinal", variable="impulse", value=nrow(impulse_data)))

##Generate bar plot
pdf(file = "deg_impulse_number_comparison.pdf", width = 7, height = 7)
p <- ggplot(data = deg_data, mapping = aes(x=pseudotime.i., y=value, fill=variable)) + 
  geom_bar(stat = "identity", position = "dodge", color="black")
p <- p + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
p <- p + xlab("Pseudotime") + ylab("Number of genes") + scale_x_discrete(labels = c("0 vs 1", "0 vs 2", "0 vs 3", "0 vs 4", "0 vs 5", "0 vs 6", "0 vs 7", "Longitudinal"))
p <- p + theme_bw() + theme(axis.text=element_text(size=14, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=12), 
                            legend.title = element_text(size=0), plot.title = element_text(size=20, face="bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
p
dev.off()
