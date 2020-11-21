library(topGO)
library(org.Hs.eg.db)
library(GO.db)
library(ComplexHeatmap)
library(RColorBrewer)
library(plyr)
library(tibble)
library(dplyr)
library(ggplot2)


#Correlation analysis between transcription and methylation--------

#Read and filter methylation data for sites
methylation_data <- read.csv("Analysis/reports_preprocessing/Rosenstiel_Control_COVID_beta_values_hg38.txt", header = TRUE, 
                             sep = '\t', row.names = 1)
meth_sample_info <- read.csv("COVID_Control_EPIC_annotation_file.csv", header = TRUE, sep = ',')

methylation_data <- methylation_data[, paste('X', as.character(meth_sample_info$Sample_ID), sep = '')]
colnames(methylation_data) <- meth_sample_info$Sample_Name2

#Read and filter transcriptome data for genes
transcriptome_data <- read.csv("../RNAseq/count_files/COVID_control_normalized_counts_070720.txt", 
                               header = TRUE, sep = '\t', row.names = 1)
transcriptome_info_file <- read.csv("../info_files/Sequenced_sample_info_020720.csv")
transcriptome_info_file <- subset(transcriptome_info_file, transcriptome_info_file$RNAseqID != '')
transcriptome_data <- transcriptome_data[, as.character(transcriptome_info_file$RNAseqID)]
colnames(transcriptome_data) <- transcriptome_info_file$Name

common_samples <- intersect(colnames(transcriptome_data), colnames(methylation_data))
transcriptome_data <- transcriptome_data[, common_samples]
methylation_data <- methylation_data[, common_samples]

#Read list of differentially expressed genes 
deg_gene_list <- read.table("../RNAseq/output/control+impulse+severe_degs.txt")
deg_gene_list <- deg_gene_list$V1

#Read and filter gene and linked methylation sites info, remove sites on sex chromosomes
gene_meth_sites <- read.csv("Analysis/Integration_with_transcriptome/gene_meth_sites_5000bp.txt", header = TRUE, 
                            sep = '\t')
gene_meth_sites <- subset(gene_meth_sites, gene_meth_sites$Chr != 'chrX' & gene_meth_sites$Chr != 'chrY' & 
                            gene_meth_sites$Chr != 'chrM')
rownames(gene_meth_sites) <- gene_meth_sites$Gene_id

deg_gene_list <- intersect(deg_gene_list, rownames(gene_meth_sites))

deg_gene_meth_sites <- gene_meth_sites[as.character(deg_gene_list), ]
deg_gene_meth_sites <- subset(deg_gene_meth_sites, deg_gene_meth_sites$no_of_meth_sites > 0)

#Calculate correlation coefficient between the gene expression of each gene and the methylation intensity of each of its link methylated sites
#Calculate FDR using a permutation approach
deg_gene_meth_sites$meth_sites <- as.character(deg_gene_meth_sites$meth_sites)
gene_meth_site_correlation <- matrix(nrow = sum(deg_gene_meth_sites$no_of_meth_sites), ncol = 5)
rownum = 1
for (i in 1:nrow(deg_gene_meth_sites)) {
  meth_sites <- strsplit(deg_gene_meth_sites[i,8], split = ';')
  #print(meth_sites)
  for (j in 1:length(meth_sites[[1]])) {
    site_id_info <- strsplit(meth_sites[[1]][j], split = ':')
    site_id <- site_id_info[[1]][1]
    distance <- site_id_info[[1]][3]
    #print(site_id)
    corr_data_frame <- data.frame(meth=t(methylation_data[site_id,]), 
                                  expr=t(transcriptome_data[as.character(deg_gene_meth_sites[i,1]),]))
    corr_data_frame <- corr_data_frame[complete.cases(corr_data_frame),]
    rho <- cor(corr_data_frame[,1], corr_data_frame[,2], method="spearman")
    fdr <- 0
    for (k in 1:1000) {
      x <- sample(corr_data_frame[,1])
      y <- sample(corr_data_frame[,2])
      rand_rho <- cor(x, y, method="spearman")
      if(abs(rand_rho) >= abs(rho)){
        fdr <- fdr + 1
      }
    }
    fdr <- fdr/1000
    gene_meth_site_correlation[rownum, ] <- c(as.character(deg_gene_meth_sites[i,1]), site_id, distance, rho, fdr)
    rownum <- rownum + 1
  }
}
gene_meth_site_correlation <- as.data.frame(gene_meth_site_correlation)
colnames(gene_meth_site_correlation) <- c("Gene", "Site", "Distance_from_TSS", "Rho", "FDR")
gene_meth_site_correlation$FDR <- as.numeric(as.character(gene_meth_site_correlation$FDR))

write.table(gene_meth_site_correlation, file = "Analysis/Integration_with_transcriptome/DEG_meth_site_correlation_5000bp_1000rep.txt", quote = FALSE,
            sep = '\t')

#Filter genes and sites for DEGs and DMPs
control_degs <- read.table("../RNAseq/output/DESeq2_results/Control_vs_COVID/control_degs.txt")$V1
longitudinal_degs <- read.table("../RNAseq/output/Longitudinal_degs.txt")$V1
control_dmps <- read.table("Analysis/Control_vs_COVID_analysis/control_dmps.txt")$V1
longitudinal_dmps <- read.table("Analysis/COVID_longitudinal_analysis/longitudinal_dmps.txt")$V1

deg_dmp_site_correlation <- subset(gene_meth_site_correlation, gene_meth_site_correlation$Gene %in% c(control_degs, longitudinal_degs)
                                           & gene_meth_site_correlation$Site %in% c(control_dmps, longitudinal_dmps))

deg_dmp_site_correlation <- unique(deg_dmp_site_correlation)
deg_dmp_site_correlation <- deg_dmp_site_correlation[order(deg_dmp_site_correlation$FDR), ]
write.table(deg_dmp_site_correlation, "Analysis/Integration_with_transcriptome/DEG_DMP_site_correlation.txt", quote = FALSE, sep = '\t', row.names = FALSE)


#GO term enrichment analysis of DMP-DEGs--------
#Extract significant DMP-DEGs
all_genes <- as.character(unique(deg_dmp_correlation$Gene))
deg_dmp_correlation <- subset(deg_dmp_correlation, deg_dmp_correlation$FDR < 0.05)
corr_genes <- as.character(unique(deg_dmp_correlation$Gene))

#Identify genes upregulated or downregulated in most pseudotimes compared to controls
pseudotime_lfc <- read.csv("../RNAseq/output/DESeq2_results/Control_vs_COVID/all_genes_pseudotime_lfc.txt", sep = '\t')
pseudotime_lfc$Gene_ID <- as.character(pseudotime_lfc$Gene_ID)
corr_genes_pseudotime_lfc <- subset(pseudotime_lfc, pseudotime_lfc$Gene_ID %in% corr_genes)
upgenes <- c()
downgenes <- c()
for (gene in corr_genes) {
  gene_pseudotime_lfc <- corr_genes_pseudotime_lfc[corr_genes_pseudotime_lfc$Gene_ID ==gene, ]
  up_count <- sum(gene_pseudotime_lfc$direction == 'up')
  down_count <- sum(gene_pseudotime_lfc$direction == 'down')
  if (up_count > down_count) {
    upgenes <- union(upgenes, gene)
  } else {
    downgenes <- union(downgenes, gene)
  }
}

#Perform GO term analysis for upregulated and downregulated genes
direction <- c("up", "down")
gene_set <- c(upgenes, downgenes)

for (i in 1:length(direction)) {
  #Define background set and test set
  inUniverse <- all_genes %in% all_genes
  inSelection <- all_genes %in% gene_set[i]
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- all_genes[inUniverse]
  
  #Run topGO analysis
  tgd <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  go_results <- GenTable(tgd, Fisher.elim = resultTopGO.elim, 
                         Fisher.classic = resultTopGO.classic,
                         orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))
  
  #Extract significant results and add genes contributing to each GO term
  go_results <- subset(go_results, go_results$Fisher.elim < 0.05)
  for(i in 1:nrow(go_results)){
    go_id <- as.vector(go_results[i,1])
    alleges <- get(go_id, org.Hs.egGO2ALLEGS)
    genes <- unlist(mget(alleges, org.Hs.egENSEMBL))
    genes_in_cat <- intersect(as.vector(genes), as.vector(upgenes))
    gene_sym_in_cat <- as.vector(unlist(mget(unlist(mget(genes_in_cat, org.Hs.egENSEMBL2EG)), org.Hs.egSYMBOL)))
    gene_sym_in_cat_str <- ""
    
    if(length(genes_in_cat) > 0){
      for(j in 1:length(gene_sym_in_cat)){
        gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j], sep = ',')
      }
    }
    if (gene_sym_in_cat_str == ""){
      gene_sym_in_cat_str <- NA
    }
    go_results$Genes[i] <- gene_sym_in_cat_str
  }
  write.table(go_results, paste("Analysis/Integration_with_transcriptome/Correlated_", direction[i], "_DEGs_GO.txt", sep = ''), sep = '\t', quote = FALSE)
  
}


#Plot selected GO terms enriched in DMP-DEGs as in figure 4I
selected_GO <- read.csv("Selected_Correlated_up_and_down_GO.csv")

selected_GO$p <- -1*log10(selected_GO$Fisher.elim)
selected_GO$n <- selected_GO$Significant/selected_GO$Annotated
selected_GO$Term <- Term(as.vector(selected_GO$GO.ID))
selected_GO$Term2 <- reorder(selected_GO$Term, 1:nrow(selected_GO))

pdf("Selected_Correlated_up_and_down_GO.pdf", width = 7.2, height = 6)
p <- ggplot(selected_GO, mapping = aes(x=Direction, y=Term2, size=n, color=p))
p <- p + geom_point()
p <- p + scale_x_discrete(limits=c("Upregulated", "Downregulated"))
p <- p + xlab("") + ylab("GO Term")
p <- p + scale_colour_gradient(high="#990000", low="#FF9999")
p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                            axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1), 
                            axis.title = element_text(size = 20))
p
dev.off()

#Analysis for comparing DMP-DEGs to coexpression modules-----

#Read gene list of each coexpression module and combine it with dmp-deg correlation data
gene_module <- read.csv("../RNAseq/output/coexpression_analysis/Without_control_and_recovery/gene_module.txt", sep = '\t')

deg_module_dmp_correlation <- inner_join(gene_module, deg_dmp_correlation, by="Gene")
deg_module_dmp_correlation$size <- -1*log10(deg_module_dmp_correlation$FDR+0.0001)

gene_list <- unique(deg_module_dmp_correlation$Gene)
gene_correlation_list <- data.frame()

#For each DEG extract the DMP with the highest correlation coefficient
for (gene in gene_list) {
  gene_correlation <- subset(deg_module_dmp_correlation, deg_module_dmp_correlation$Gene == gene)
  gene_correlation <- gene_correlation[which(abs(gene_correlation$Rho) == max(abs(gene_correlation$Rho)))[1], ]
  gene_correlation_list <- rbind(gene_correlation_list, data.frame(gene_correlation))
}
colnames(gene_correlation_list) <- colnames(deg_module_dmp_correlation)

#Add module numbers to the table and remove the grey module
module_colors_numbers <- read.csv("../RNAseq/output/coexpression_analysis/Without_control_and_recovery/Module_name_to_number.csv", row.names = 1)
gene_correlation_list$Module <- module_colors_numbers[as.character(gene_correlation_list$ModuleColor), "Module_number"]
gene_correlation_list <- subset(gene_correlation_list, gene_correlation_list != "grey")

#Plot the correlation coefficients for DMP-DEGs in modules as in figure S5C
pdf("Analysis/Integration_with_transcriptome/Coexpression_module_DMP_correlation.pdf", width = 8, height = 6)
p <- ggplot(gene_correlation_list, aes(x=Module, y=Rho))
p <- p + geom_point(aes(size=size), position=position_jitter(width=0.2, height = 0.1), pch=21, fill="darkgrey")
p <- p + scale_x_discrete(limits=c(as.character(module_colors_numbers$Module_number)))
p <- p + ylab("Spearman's Rho") + xlab("Module")
p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), axis.text.x = element_text(size=14, color = "black", angle=45, hjust=1), 
                            axis.title = element_text(size = 20), legend.position = "bottom")
p
dev.off()

#Chi-square test to identify modules containing more or less than expected number of DMP-DEGs
significant_genes <- as.character(unique(deg_dmp_correlation$Gene))

module_names <- unique(gene_module$ModuleColor)
module_ecg_chisq <- data.frame()

for (module in module_names) {
  gene_module$in_module <- ifelse(gene_module$ModuleColor == module, "in", "out")
  gene_module$ecg <- ifelse(gene_module$Gene %in% significant_genes, "ECG", "not-ECG")
  chisq <- chisq.test(gene_module$in_module, gene_module$ecg, correct = FALSE)
  module_ecg_chisq <- rbind(module_ecg_chisq, data.frame(module, chisq$observed[1]/chisq$observed[3], chisq$observed[2]/chisq$observed[4], chisq$observed[1], 
                                                         chisq$expected[1],chisq$statistic, chisq$p.value)) 
}

colnames(module_ecg_chisq) <- c("Module", "Fraction_ECG_in_module", "Fraction_ECG_out_module", "Observed", "Expected","Chisq", "P")
module_ecg_chisq$padj <- p.adjust(module_ecg_chisq$P, method = "BH")
module_ecg_chisq$Ratio <- module_ecg_chisq$Observed/module_ecg_chisq$Expected
module_ecg_chisq <- subset(module_ecg_chisq, module_ecg_chisq$Module != "grey")
module_ecg_chisq$Modulename <- module_colors_numbers[as.character(module_ecg_chisq$Module), "Module_number"]

#Plot Observed-Expected ratio as in figure S5D
pdf("Analysis/Integration_with_transcriptome/Coexpression_module_ECG_comparison.pdf", width = 6, height = 5)
p <- ggplot(data = module_ecg_chisq, mapping = aes(x=Modulename, y=Ratio-1, label=paste("P =", signif(padj,3)))) + 
  geom_bar(stat = "identity", color="black", fill="darkgrey")
p <- p + geom_text(data=subset(module_ecg_chisq, module_ecg_chisq$Ratio > 1), size=3, nudge_y = 0.05)
p <- p + geom_text(data=subset(module_ecg_chisq, module_ecg_chisq$Ratio < 1), size=3, nudge_y = -0.05)
#p <- p + scale_fill_manual(values=levels(module_ecg_chisq$Module))
p <- p + ylab("Observed/Expected - 1")  + xlab("Module")
p <- p + scale_x_discrete(limits=as.character(module_colors_numbers$Module_number))
p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), axis.text.x = element_text(size=14, color = "black", angle=45, hjust=1), 
                            axis.title = element_text(size = 20), legend.position = "bottom")
p
dev.off()

write.table(module_ecg_chisq, "Analysis/Integration_with_transcriptome/Coexpression_module_ECG_comparison.txt", quote = FALSE, sep = '\t', row.names = FALSE)

#Cell type expression of DMP-DEGs-----
#Add gene names to the correlation table
mart_export <- read.table("../reference/human/gencode_v25_gene_info.txt", header = T, sep="\t")
unique_mart <- subset(mart_export, duplicated(mart_export$Gene_id) == FALSE)
rownames(unique_mart) <- unique_mart$Gene_id
deg_dmp_correlation$Gene_name <- unique_mart[as.character(deg_dmp_correlation$Gene), ]$Gene_name

#Filter for cannonical DMP-DEG pairs
deg_dmp_correlation <- subset(deg_dmp_correlation, deg_dmp_correlation$Rho < 0)

corr_genes <- as.character(unique(deg_dmp_correlation$Gene_name))

#Identify DEGs upregulated in the inflammatory pseusotimes (early) or downregulated in the convalescent pseusotimes (late)
pseudotime_lfc <- read.csv("../RNAseq/output/DESeq2_results/Control_vs_COVID/all_genes_pseudotime_lfc.txt", sep = '\t')
pseudotime_lfc$Gene_ID <- as.character(pseudotime_lfc$Gene_ID)
corr_genes_pseudotime_lfc <- subset(pseudotime_lfc, pseudotime_lfc$gene %in% corr_genes)
upgenes_early <- c()
downgenes_late <- c()
for (gene in corr_genes) {
  gene_pseudotime_lfc <- corr_genes_pseudotime_lfc[corr_genes_pseudotime_lfc$gene ==gene, ]
  up_count <- sum(gene_pseudotime_lfc[gene_pseudotime_lfc$pseudotime %in% c(1,2,3,4), ]$direction == 'up')
  down_count <- sum(gene_pseudotime_lfc[gene_pseudotime_lfc$pseudotime %in% c(5,6), ]$direction == 'down')
  if (up_count > 0) {
    upgenes_early <- union(upgenes_early, gene)
  } 
  if (down_count > 0) {
    downgenes_late <- union(downgenes_late, gene)
  }
}

#Add coexpression module to the table
gene_module$Module <- module_color_to_number[as.character(gene_module$ModuleColor), "Module_number"]
gene_module$ModuleLabel <- NULL
combined_data <- left_join(deg_dmp_correlation, gene_module, by=c("Gene"="Gene"))
combined_data$Module <- as.character(combined_data$Module)
combined_data[is.na(combined_data$Module), "Module"] <- "none"

#Cell type expression analysis for genes upregulated in the inflammatory pseudotimes
gene_list <- upgenes_early

#Get sample list for scRNAseq samples
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
matDataHeat_gene_module <- combined_data[match(rownames(matDataHeat), combined_data$Gene_name.x), c("Gene_name.x", "Module", "ModuleColor")]
matDataHeat_gene_module[is.na(matDataHeat_gene_module$ModuleColor), "ModuleColor"] <- "grey"
matDataHeat_gene_module$Module <- factor(matDataHeat_gene_module$Module, levels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9","M10","none"))
matDataHeat_gene_module <- matDataHeat_gene_module[order(matDataHeat_gene_module$Module), ]
matDataHeat <- matDataHeat[as.character(matDataHeat_gene_module$Gene_name), ]

#Plot heatmap as in figure S5E
col = as.character(unique(matDataHeat_gene_module$ModuleColor))
names(col) <- unique(matDataHeat_gene_module$Module)


row_ha <- rowAnnotation(df = data.frame(Module = matDataHeat_gene_module$Module),
                        col = list(Module = col))


pdf("Analysis/Integration_with_transcriptome/Negatively_correlated_DEG_up_in_early_module_cell_type_heatmap.pdf", width = 8, height = 16)
Heatmap(matDataHeat, col=rev(brewer.pal(11, "RdYlBu")), 
        row_names_gp = gpar(fontsize = 1), left_annotation = row_ha,
        cluster_columns = FALSE, column_names_side = "top", column_names_rot = 45, 
        column_split = factor(colnames(matDataHeat), levels = colnames(matDataHeat)), 
        column_title_gp = gpar(fontsize=0),
        row_split = matDataHeat_gene_module$Module, border = TRUE, cluster_row_slices = FALSE, row_title_rot = 0)
dev.off()

#Cell type expression analysis for genes downregulated in the convalescence pseudotimes
gene_list <- downgenes_late

#Get sample list for scRNAseq samples
sample_list <- read.csv("../info_files/Sequenced_sample_info_020720.csv")
sample_list <- subset(sample_list, sample_list$Pseudotime %in% c(5,6))
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
matDataHeat_gene_module <- combined_data[match(rownames(matDataHeat), combined_data$Gene_name.x), c("Gene_name.x", "Module", "ModuleColor")]
matDataHeat_gene_module[is.na(matDataHeat_gene_module$ModuleColor), "ModuleColor"] <- "grey"
matDataHeat_gene_module$Module <- factor(matDataHeat_gene_module$Module, levels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9","M10","none"))
matDataHeat_gene_module <- matDataHeat_gene_module[order(matDataHeat_gene_module$Module), ]
matDataHeat <- matDataHeat[as.character(matDataHeat_gene_module$Gene_name), ]

#Plot heatmap as in figure S5F
col = as.character(unique(matDataHeat_gene_module$ModuleColor))
names(col) <- unique(matDataHeat_gene_module$Module)


row_ha <- rowAnnotation(df = data.frame(Module = matDataHeat_gene_module$Module),
                        col = list(Module = col))

pdf("Analysis/Integration_with_transcriptome/Negatively_correlated_DEG_down_in_late_module_cell_type_heatmap.pdf", width = 8, height = 10)
Heatmap(matDataHeat, col=rev(brewer.pal(11, "RdYlBu")), 
        row_names_gp = gpar(fontsize = 1), left_annotation = row_ha,
        cluster_columns = FALSE, column_names_side = "top", column_names_rot = 45, 
        column_split = factor(colnames(matDataHeat), levels = colnames(matDataHeat)), 
        column_title_gp = gpar(fontsize=0),
        row_split = matDataHeat_gene_module$Module, border = TRUE, cluster_row_slices = FALSE, row_title_rot = 0)
dev.off()

