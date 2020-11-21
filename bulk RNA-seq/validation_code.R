library(DESeq2)
library(WGCNA)
library(reshape2)
library(tidyverse)
library(ggplot2)


#Code to compare module eigengenes between sampling points in cohort 3------

#Read gene list for cohort 1 modules
gene_list <- read.csv("../RNAseq/output/coexpression_analysis/Without_control_and_recovery/gene_module.txt", sep = '\t')

#Read sample list for cohort 3
sample_list <- read.csv("info_files/validation_samples_info_250720.csv")
sample_list$RUMC.ID <- as.character(sample_list$RUMC.ID)

#Read count data for cohort 3
count_data <- read.csv("count_files/raw_counts2020-07-10.txt", sep = '\t')
rownames(count_data) <- gsub('(.*)\\..*', '\\1', rownames(count_data))

#Edit sample names and order count data 
sample_list$sample_names <- paste('X', sample_list$RUMC.ID, sep = '')
sample_list$predeath_sample_name <- paste(sample_list$PatientID, sample_list$Death, sample_list$Predeath_timepoints, sep = '_')
count_data <- count_data[, as.character(sample_list$sample_names)]

#Prepare col data
col_data <- sample_list[, c("PatientID", "sample_names", "Death", "Age", "Predeath_timepoints", 
                            "predeath_sample_name")]
col_data$Death <- as.factor(as.character(col_data$Death))

#Create DESeq object and normalize
dds <- DESeqDataSetFromMatrix(count_data, colData = col_data, design = ~ PatientID + Predeath_timepoints)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
normalized_expression_counts <- counts(dds, normalized=TRUE)

#Get common genes between both cohorts
gene_list <- subset(gene_list, gene_list$Gene %in% rownames(normalized_expression_counts))

#Calculate module eigengenes from cohort 3 data
sig_genes_expression_data <- normalized_expression_counts[as.character(gene_list$Gene),]
log_counts <- log2(sig_genes_expression_data + 1)
gene_list$ModuleColor <- as.character(gene_list$ModuleColor)
MEsO <- moduleEigengenes(t(log_counts), gene_list$ModuleColor)$eigengenes
MEs <- orderMEs(MEsO)
rownames(MEs) <- rownames(t(log_counts))

#Add metadata to module eigengenes
MEs[, c("PatientID", "Death", "Predeath_timepoints", "predeath_sample_name")] <- 
  sample_list[, c("PatientID", "Death", "Predeath_timepoints", "predeath_sample_name")]
MEs$Death <- ifelse(MEs$Death == 0, "Survived", "Died")

#Perform paired Mann Whitney U test to compare eigengenes between sampling points in survivors and non-survivors
melted_ME <- melt(MEs, id.vars = c("PatientID", "Death", "Predeath_timepoints", "predeath_sample_name"))
wilcox_summary <- melted_ME %>% 
  group_by(variable, Death) %>%
  summarise(p=round(wilcox.test(value ~ Predeath_timepoints, paired=TRUE)$p.value,2),
            y=median(value),
            x=1) 
wilcox_summary <- as.data.frame(wilcox_summary)

write.table(wilcox_summary, "output/module_preservation/Death_ME_comparison_wilcox.txt", sep = '\t', quote = FALSE, row.names = FALSE)

#Plot module eigengene values for M2, M4 and M7 as in figure 7B

plot_data <- melt(MEs, id.vars = c("Death", "Predeath_timepoints"), measure.vars = c("MEfloralwhite", "MEplum1", "MEbrown"))
plot_data$Module <- ifelse(plot_data$variable=="MEfloralwhite", "M2", ifelse(plot_data$variable=="MEplum1", "M4", "M7"))
plot_data$x <- ifelse(plot_data$Death == "Died", "Non-survivor", "Survivor")
plot_data$col <- paste(plot_data$x, plot_data$Predeath_timepoints, sep = '_')

cairo_ps("output/module_preservation/Death_ME_comparison.ps", width=10, height = 4)
p <- ggplot(data = plot_data, aes(y=value, x=x, fill=col)) + geom_boxplot()
p <- p + facet_wrap(~Module, nrow=1, scale="free")
p <- p + ylab("Module Eigengene") + xlab("")
p <- p + scale_fill_manual(values = c("#FF6666", "#990000", "#66B2FF", "#000099"))
p <- p + theme_minimal() + theme(axis.text=element_text(size=14, color = "black"),axis.title=element_text(size=16), 
                                 axis.line.y = element_line(size=0.5), strip.text.x = element_text(size=14, face = "bold"),
                                 legend.text = element_text(size=12), legend.title = element_text(size=0), plot.title = element_text(size=20, face="bold", hjust = 0.5))
p
dev.off()



#Code for differential expression analysis between sampling points in cohort 3------
#Read gene info file
mart_export <- read.table("../reference/human/gencode_v25_gene_info.txt", header = T, sep="\t")
unique_mart <- subset(mart_export, duplicated(mart_export$Gene_id) == FALSE)
rownames(unique_mart) <- unique_mart$Gene_id

#Classify survivors and non survivors
death_sample_list <- subset(sample_list, sample_list$Death == 1)
no_death_sample_list <- subset(sample_list, sample_list$Death == 0)

#Deseq analysis for non survivor samples
death_count_data <- count_data[, as.character(death_sample_list$sample_names)]
col_data <- death_sample_list[, c("PatientID", "sample_names", "Death","Age", "Predeath_timepoints", 
                                  "predeath_sample_name", "predeath_sample_name")]
col_data$Death <- as.factor(as.character(col_data$Death))
dds_count <- DESeqDataSetFromMatrix(death_count_data, colData = col_data, design = ~ PatientID + Predeath_timepoints)
dds_count <- dds_count[ rowSums(counts(dds_count)) > 1, ]
dds_count <- estimateSizeFactors(dds_count)
dds <- DESeq(dds_count, betaPrior = FALSE)
res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
res_sorted <- res[order(res$padj), ]
output_folder <- "output/DESeq2_results/Death_A_vs_B/"
death_res_sorted <- res_sorted
death_res_sorted$gene <- unique_mart[rownames(death_res_sorted), ]$Gene_name
write.table(death_res_sorted, file = file.path(output_folder, "DESeq2result_genenames_A_vs_B.txt"), 
            sep="\t", quote=FALSE)

#Deseq analysis for survivor samples
no_death_count_data <- count_data[, as.character(no_death_sample_list$sample_names)]
col_data <- no_death_sample_list[, c("PatientID", "sample_names", "Death", "Age", "Predeath_timepoints", 
                                     "predeath_sample_name", "predeath_sample_name")]
col_data$Death <- as.factor(as.character(col_data$Death))
dds_count <- DESeqDataSetFromMatrix(no_death_count_data, colData = col_data, design = ~ PatientID + Predeath_timepoints)
dds_count <- dds_count[ rowSums(counts(dds_count)) > 1, ]
dds_count <- estimateSizeFactors(dds_count)
dds <- DESeq(dds_count, betaPrior = FALSE)
res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
res_sorted <- res[order(res$padj), ]
output_folder <- "output/DESeq2_results/no_Death_A_vs_B/"
no_death_res_sorted <- res_sorted
no_death_res_sorted$gene <- unique_mart[rownames(no_death_res_sorted), ]$Gene_name
write.table(no_death_res_sorted, file = file.path(output_folder, "DESeq2result_genenames_A_vs_B.txt"), 
            sep="\t", quote=FALSE)

#Code for comparison between cohort 1 DEGs and cohort 3 non survivor DEGs------
#Compare DEGs in non survivors to cohort 1 gene list
module_color_to_number <- read.csv("../RNAseq/output/coexpression_analysis/Without_control_and_recovery/Module_name_to_number.csv", row.names = 1)
gene_module$Module <- module_color_to_number[as.character(gene_module$ModuleColor), "Module_number"]
gene_module$ModuleLabel <- NULL

death_res_sorted$Gene_ID <- rownames(death_res_sorted)
combined_data <- left_join(death_res_sorted, gene_module, by=c("Gene_ID"="Gene"))
combined_data$log_padj <- -1 * log10(combined_data$padj)
combined_data$significant <- ifelse(combined_data$padj < 0.05, "yes", "no")
combined_data$significant <- as.factor(combined_data$significant)

#Plot volcano plot for deseq result from non survivors and color genes by cohort 1 modules as in figure 7C
cairo_ps("output/DESeq2_results/Death_A_vs_B/Volcano_plot_colored_by_module.ps")
p <- ggplot(data = combined_data, aes(x=log2FoldChange, y=log_padj, fill=Module, alpha=significant)) + 
  geom_point(pch=21, size=3, color="black") 
p <- p + geom_vline(xintercept = 0, lty=3) + geom_vline(xintercept = 1, lty=4) + geom_vline(xintercept = -1, lty=4)
p <- p + geom_hline(yintercept = 1.30103, lty=4) + geom_hline(yintercept = 3, lty=4)
p <- p + scale_fill_manual(values=rownames(module_color_to_number), 
                           limits = as.character(module_color_to_number$Module_number), na.value="darkgrey")
p <- p + scale_alpha_manual(values = c(0.2, 1))
p <- p + xlab("Expression (Log2 Fold Change)") + ylab("- Log10 Q-value")
p <- p + theme_bw() + theme(axis.text = element_text(size = 14, color="black"), axis.title = element_text(size = 16))
p
dev.off()

#Fisher exact test to determine the significance of overlap between non survivor DEGs and each of the cohort 1 modules
module_names <- unique(module_color_to_number$Module_number)
module_sig_fisher <- data.frame()
combined_data$Module <- as.character(combined_data$Module)
combined_data[is.na(combined_data$Module), "Module"] <- "none"
for (module in module_names) {
  in_module <- ifelse(combined_data$Module == module, "in", "out")
  significant <- combined_data$significant
  significant <- factor(significant, levels=c("yes", "no"))
  fisher <- fisher.test(table(in_module, significant))
  module_sig_fisher <- rbind(module_sig_fisher, data.frame(module, fisher$estimate, fisher$p.value))
}

colnames(module_sig_fisher) <- c("Module", "Odds_ratio", "P")
module_sig_fisher$padj <- p.adjust(module_sig_fisher$P, method = "BH")

write.table(module_sig_fisher, "output/module_preservation/module_significant_in_death_fisher_results.txt", sep = '\t', quote = FALSE, row.names = FALSE)


#Cell type expression of cohort 3 non survivor DEGs in cohort 1-------
#Extract non survivor DEGs
combined_data <- subset(combined_data, combined_data$padj < 0.05)
combined_data$Module <- as.character(combined_data$Module)
combined_data[is.na(combined_data$Module), "Module"] <- "none"
gene_list <- as.character(combined_data$gene)
gene_list <- gene_list[!is.na(gene_list)]

#Get sample list for scRNAseq samples and filter for pseudotimes 1, 2 and 3
sample_list <- read.csv("../info_files/Sequenced_sample_info_020720.csv")
sample_list <- subset(sample_list, sample_list$Pseudotime %in% c(1,2,3))
sample_list <- subset(sample_list, sample_list$scRNAseqID != '')
sample_list$sample_name <- paste(sample_list$Name, sample_list$Pseudotime, sep = '_')

#Read average cell type data and filter for gene_list
cell_types <- c("Monocytes", "Granulocytes", "Erythroid", "NK", "Proliferative", "CD4", "CD8", "DC", "Bcells", "Plasmablasts", "Precursors", 
                "Megakaryocytes")
sample_names <- read.csv("../RNAseq/output/Intergration_with_sc_data/Merge_sample_average.csv", row.names = 2)
folder <- "../RNAseq/output/Intergration_with_sc_data/"
all_cell_type_data <- data.frame()
for (cell in cell_types) {
  cell_type_data <- read.csv(file.path(folder, paste(cell, "average.csv", sep = '_')), row.names = 1)
  colnames(cell_type_data) <- as.character(sample_names[colnames(cell_type_data), ])
  cell_type_data <- cell_type_data[intersect(gene_list, rownames(cell_type_data)), 
                                   intersect(as.character(sample_list$Name), colnames(cell_type_data))]
  colnames(cell_type_data) <- sample_list[match(colnames(cell_type_data), sample_list$Name), "sample_name"]
  colnames(cell_type_data) <- paste(colnames(cell_type_data), cell, sep = '_')
  all_cell_type_data <- merge(all_cell_type_data, cell_type_data, by=0, all=TRUE)
  all_cell_type_data <- column_to_rownames(all_cell_type_data, var="Row.names")
}

#Avergae expression for each across pseudotimes cell type
ct <- str_split_fixed(colnames(all_cell_type_data), '_', 5)[, 5]
matDataHeat <- do.call(cbind, lapply(cell_types, function(tp) {
  vecidxCols <- which(ct %in% tp)
  if (length(vecidxCols) > 1) {
    return(rowMeans(all_cell_type_data[, vecidxCols], na.rm = TRUE))
  } else {
    return(all_cell_type_data[, vecidxCols])
  }
}))
colnames(matDataHeat) <- cell_types

#Remove genes with undetectable expression and scale data
matDataHeat <- matDataHeat [which(rowSums(matDataHeat) > 0), ]
matDataHeat <- t(scale(t(matDataHeat)))

#Assign modules to heatmap genes
matDataHeat_gene_module <- combined_data[match(rownames(matDataHeat), combined_data$gene), c("gene", "Module", "ModuleColor")]
matDataHeat_gene_module[is.na(matDataHeat_gene_module$ModuleColor), "ModuleColor"] <- "grey"
matDataHeat_gene_module$Module <- factor(matDataHeat_gene_module$Module, levels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M10","none"))
matDataHeat_gene_module <- matDataHeat_gene_module[order(matDataHeat_gene_module$Module), ]
matDataHeat <- matDataHeat[as.character(matDataHeat_gene_module$gene), ]

#Plot heatmap as in figure 7D
col = as.character(unique(matDataHeat_gene_module$ModuleColor))
names(col) <- unique(matDataHeat_gene_module$Module)


row_ha <- rowAnnotation(df = data.frame(Module = matDataHeat_gene_module$Module),
                        col = list(Module = col), 
                        gp = gpar(col = "black"))

pdf("output/module_preservation/Death_deg_module_cell_type_heatmap.pdf", width = 8, height = 14)
Heatmap(matDataHeat, col=rev(brewer.pal(11, "RdYlBu")), 
        row_names_gp = gpar(fontsize = 6), left_annotation = row_ha,
        cluster_columns = FALSE, column_names_side = "top", column_names_rot = 45, 
        column_split = factor(colnames(matDataHeat), levels = colnames(matDataHeat)), 
        column_title_gp = gpar(fontsize=0),
        row_split = matDataHeat_gene_module$Module, border = TRUE, cluster_row_slices = FALSE, row_title_rot = 0)
dev.off()
