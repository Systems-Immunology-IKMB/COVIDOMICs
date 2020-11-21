library(DESeq2)
library(ggplot2)
library(magrittr)
library(WGCNA)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(GO.db)
library(dplyr)
library(LOLA)
library(GenomicRanges)

setwd("/Users/neha/Projects/ProjectC/RNAseq/")

folder <- "output/coexpression_analysis/Without_control_and_recovery//"

if (!dir.exists(folder)){
  dir.create(folder)
}

exclude_pseudotimes <- c('2', '7')

#Code for coexpression modules generation ---------------------------
# Preprocess data
covid_info_file <- read.csv("../info_files/Sequenced_sample_info_020720.csv")
covid_info_file <- subset(covid_info_file, covid_info_file$RNAseqID != '')

covid_info_file$Sample_name <- paste(covid_info_file$Name, covid_info_file$Pseudotime, sep = '_')

normalized_expression_measures <- read.csv("count_files/COVID_control_normalized_counts_070720.txt", sep = '\t')
normalized_expression_measures <- normalized_expression_measures[, as.character(covid_info_file$RNAseqID)]
colnames(normalized_expression_measures) <- covid_info_file$Sample_name
sig_genes <- read.csv("output/control+impulse+severe_degs.txt", header = FALSE)
sig_genes <- sig_genes$V1

sig_genes_expression_data <- normalized_expression_measures[as.character(sig_genes),]

log_counts <- log2(sig_genes_expression_data + 1)

#Filter for required pseudotimes
filtered_covid_info_file <- subset(covid_info_file, !(covid_info_file$Pseudotime %in% exclude_pseudotimes))

#Define wgcna_matrix
wgcna_matrix <- log_counts[, as.character(filtered_covid_info_file$Sample_name)]

wgcna_matrix <- t(wgcna_matrix)
s <- bicor(wgcna_matrix)

#Select parameters for scale free topology
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(t(log_counts), powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

#Create adjacency matrix
adj_mat <- adjacency.fromSimilarity(s, power=14, type = 'signed')
TOM = TOMsimilarity(adj_mat)
dissTOM = 1-TOM

#Make gene tree
gene_tree <- hclust(as.dist(1 - adj_mat), method = "average")
TOM_gene_tree <- hclust(as.dist(dissTOM), method = "average")

#Define modules

module_labels <- cutreeDynamicTree(dendro=TOM_gene_tree, minModuleSize=15, 
                                   deepSplit=TRUE)
module_colors <- labels2colors(module_labels)

#Merge modules
mergedColor <- mergeCloseModules(wgcna_matrix,module_colors,cutHeight=.45)$color


#Output modules
gene_module <- data.frame(Gene = colnames(wgcna_matrix), ModuleLabel = module_labels, ModuleColor = mergedColor)
mart_export <- read.csv("../..//reference/human/gencode_v25_gene_info.txt", header = T, sep="\t")
unique_mart <- subset(mart_export, duplicated(mart_export$Gene_id) == FALSE)
rownames(unique_mart) <- unique_mart$Gene_id
rownames(gene_module) <- gene_module$Gene
gene_module$Gene_name <- unique_mart[rownames(gene_module), ]$Gene_name
write.table(gene_module, file.path(folder, "gene_module.txt"), 
            quote = FALSE, sep = '\t', row.names = FALSE)

#Code to calculate module eigengenes and plot them as in figure 4A and figure S3F------------
#Calculate module eigengenes for unfiltered dataset
MEsO <- moduleEigengenes(t(log_counts), mergedColor)$eigengenes
MEs <- orderMEs(MEsO)
rownames(MEs) <- rownames(t(log_counts))

write.table(MEs, file.path(folder, "module_eigengenes.txt"), 
            quote = FALSE, sep = '\t')

eigengenes <- as.data.frame(t(MEs))
filtered_eigengenes <- eigengenes[!(rownames(eigengenes)=="MEgrey"), ]

module_colors <- gsub("ME", "", rownames(filtered_eigengenes))

#Average eigengene values for each pseudotime

vecTP <- str_split_fixed(colnames(pt_sorted_eigengenes), '_', 4)[, 4]
vecUniqueTP <- unique(vecTP)

matDataHeat <- do.call(cbind, lapply(vecUniqueTP, function(tp) {
  vecidxCols <- which(vecTP %in% tp)
  if (length(vecidxCols) > 1) {
    return(rowMeans(pt_sorted_eigengenes[, vecidxCols], na.rm = TRUE))
  } else {
    return(pt_sorted_eigengenes[, vecidxCols])
  }
}))

colnames(matDataHeat) <- vecUniqueTP
rownames(matDataHeat) <- gsub('ME','',rownames(matDataHeat))
module_color_names <- read.csv(file.path(folder, "Module_name_to_number.csv"), row.names = 1)
module_names <- module_color_names[rownames(matDataHeat), ]
rownames(matDataHeat) <- module_names

#Eigengene heatmap as in figure 4A

col_ha = HeatmapAnnotation(df = data.frame(Pseudotime = colnames(matDataHeat)),
                           col = list(Pseudotime = c('0' = '#A7A9AC', '1' = '#E84F8C', '2' = '#7031B9', '3' = '#A04E9E', 
                                                     '4' = '#E65826', '5' = '#F99B1C', '6' = '#FDC077', '7' = '#51BBFE')), 
                           gp = gpar(col = "black"))
bar_data <- table(gene_module$ModuleColor)
bar_data <- bar_data[module_colors]
names(bar_data) <- rownames(matDataHeat)

row_ha2 = rowAnnotation('#genes' = anno_barplot(as.vector(bar_data), which = "row", ylim=c(0,1800))) 
ha_barplot <- HeatmapAnnotation(barplot = anno_barplot(as.vector(bar_data), which = "row"), 
                                which = "row", width = unit(2, "cm"))

pdf(file.path(folder, "Eigengene_heatmap1.pdf"), width=6, height=5)
Heatmap(as.matrix(matDataHeat), col=rev(brewer.pal(11, "RdYlBu")), top_annotation = col_ha, 
        right_annotation = row_ha2, cluster_columns = FALSE)
dev.off()

#Eigengene heatmap as in figure S3F 
Age <- covid_info_file$Age
Sex <- covid_info_file$Gender
Pseudotime <- as.character(covid_info_file$Pseudotime)
Patient <- covid_info_file$PatientID

n <- length(unique(Patient))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = col_vector[1:n]
names(col_vector) <- as.character(unique(Patient))

sorted_info <- covid_info_file[order(covid_info_file$Pseudotime), ]
pt_sorted_eigengenes <- filtered_eigengenes[, as.character(sorted_info$Sample_name)]

module_names <- module_color_names[gsub('ME', '',rownames(pt_sorted_eigengenes)), ]
rownames(pt_sorted_eigengenes) <- as.character(module_names)

col_ha = HeatmapAnnotation(df = data.frame(age = sorted_info$Age, gender = sorted_info$Gender, 
                                           patient = sorted_info$PatientID, pseudotime = as.character(sorted_info$Pseudotime)),
                           col = list(age = colorRamp2(c(min(Age), max(Age)), c("white", "orange")),
                                      gender = c('m' = '#009999', 'f' = '#7F00FF'),
                                      patient = col_vector, 
                                      pseudotime = c('0' = '#A7A9AC', '1' = '#E84F8C', '2' = '#7031B9', '3' = '#A04E9E', 
                                                     '4' = '#E65826', '5' = '#F99B1C', '6' = '#FDC077', '7' = '#51BBFE')))

pdf(file.path(folder, "Eigengene_heatmap2.pdf"), width=13, height=5)
Heatmap(as.matrix(pt_sorted_eigengenes), col=rev(brewer.pal(11, "RdYlBu")), top_annotation = col_ha, show_column_names = FALSE, 
        cluster_columns = FALSE)
dev.off()



#Code for gene ranking based on module membership used for GSEA ---------------------------

module_colors <- unique(gene_module$ModuleColor)

#Calculate the module membership scores
MEs <- MEs[rownames(wgcna_matrix), ]
gene_module_relationship <- as.data.frame(cor(wgcna_matrix, MEs, use='p'))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(gene_module_relationship), nrow(wgcna_matrix)))

write.table(gene_module_relationship, file.path(folder, "module_membership.txt"), sep = '\t', quote = FALSE)

gene_module_relationship$Gene <- unique_mart[rownames(gene_module_relationship), ]$Gene_name
gene_module_relationship <- gene_module_relationship[, c(length(module_colors) + 1, 1:length(module_colors))]

write.table(gene_module_relationship, file.path(folder, "gene_symbol_module_membership.txt"), sep = '\t', quote = FALSE, row.names = FALSE)

#Rank genes by module membership score for each module and output
for (i in 1:length(module_colors)) {
  module_gmm <- data.frame(Gene=gene_module_relationship$Gene, Rank=gene_module_relationship[[paste('ME', module_colors[i], sep = '')]])
  module_gmm <- module_gmm[order(module_gmm[,2], decreasing = TRUE), ]
  write.table(module_gmm, file.path(folder, paste("GSEA/GMM_rank_module_", module_colors[i], ".rnk", sep = '')), sep = '\t', row.names = FALSE, quote = FALSE)
}




#Code for calculating correlation between module eigengenes and clinical parameters as well as single cell fractions as in figures 4B and 4C--------
#Read clinical data file
clinical_data <- read.csv("../info_files/Clinical_data_sampling_points.csv")
rownames(clinical_data) <- clinical_data$Name
clinical_data <- clinical_data[, c("Ordinal.scale", "TimepointFromOnset", "HR.heart.rate.", "Thrombopoetin", "D.Dimers", 
                                   "White.blood.cells", "Neutrophils", "Lymphocytes", "Hematokrit", "Platelet", "CRP", "IL6", 
                                   "Ferritin", "Glucose", "Troponin.T", "antiviral.IgA", "antiviral.IgG")]
#Remove grey module
ME <- ME[, colnames(ME)!="MEgrey"]

#Get common samples having both clinical data and module eigengenes
rownames(ME) <- gsub('_[0,1,2,3,4,5,6,7]$', '', rownames(ME))
common_samples <- intersect(rownames(ME), rownames(clinical_data))

#Filter eigengenes and clinical data for common samples
filtered_ME <- ME[common_samples, ]
filtered_clinical_data <- clinical_data[common_samples, ]

#Calculate correlation between module eigenegens and clinical data and output
module_trait_cor <- cor(filtered_ME, filtered_clinical_data, use='p', method = "spearman")
module_trait_cor_pvalue <- corPvalueStudent(module_trait_cor, nrow(filtered_ME))
module_trait_cor_adj_pvalue <- p.adjust(module_trait_cor_pvalue, method = 'BH')
module_trait_cor_adj_pvalue <- matrix(module_trait_cor_adj_pvalue, nrow = nrow(module_trait_cor_pvalue), ncol = ncol(module_trait_cor_pvalue))
colnames(module_trait_cor_adj_pvalue) <- colnames(module_trait_cor_pvalue)
rownames(module_trait_cor_adj_pvalue) <- rownames(module_trait_cor_pvalue)
write.table(module_trait_cor, 
            file = file.path(folder, "module_trait_correlation_240720.txt"), 
            sep = '\t', quote = FALSE)
write.table(module_trait_cor_adj_pvalue, 
            file = file.path(folder, "module_trait_correlation_P_240720.txt"), 
            sep = '\t', quote = FALSE)

#Edit rownames and colnames for output

rownames(module_trait_cor) <- gsub('ME', '', rownames(module_trait_cor))
rownames(module_trait_cor_adj_pvalue) <- gsub('ME', '', rownames(module_trait_cor_adj_pvalue))
module_trait_cor <- module_trait_cor[rownames(module_color_names), ]
module_trait_cor_adj_pvalue <- module_trait_cor_adj_pvalue[rownames(module_color_names), ]
rownames(module_trait_cor) <- as.character(module_color_names$Module_number)
rownames(module_trait_cor_adj_pvalue) <- as.character(module_color_names$Module_number)
colnames(module_trait_cor) <- c("Ordinal scale", "Timepoint from onset", "Heart rate", "Thrombopoetin", " D-Dimers", "White blood cells", 
                                "Neutrophils", "Lymphocytes", "Hematokrit", "Platelet", "CRP", "IL-6", "Ferritin", "Glucose", "Troponin T", 
                                "anti SARS-CoV-2 IgA", "anti SARS-CoV-2 IgG")

#Plot correlation heatmap as in figure 4C
pdf(file = file.path(folder, "module_trait_correlation_250820.pdf"), 
    width=12, height = 10)
corrplot(module_trait_cor, method="color", type="full", col=viridis::viridis(30), p.mat = module_trait_cor_adj_pvalue, sig.level = c(.001, .01, .05), pch.cex=1.5, pch.col="white", 
         tl.col = "black", insig="label_sig", tl.srt = 45)
dev.off()

#Read scRNAseq cell numbers file
cell_numbers <- read.csv("output/Intergration_with_sc_data/Cell_numbers_scRNAseq_230720.csv")
rownames(cell_numbers) <- cell_numbers$Sample
cell_numbers[, c("Sample", "Disease", "Total.cells", "Total.cells2")] <- NULL
ind <- apply(cell_numbers, 1, function(x) all(is.na(x)))
cell_numbers <- cell_numbers[!ind, ]
cell_fraction <- cell_numbers/rowSums(cell_numbers)

#Get common samples having both scRNAseq data and module eigengenes
common_samples <- intersect(rownames(ME), rownames(cell_numbers))

#Filter eigengenes and cell fraction for common samples
filtered_ME <- ME[common_samples, ]
filtered_cell_fraction <- cell_fraction[common_samples, ]

#Calculate correlation between module eigenegens and cell fractions and output
module_trait_cor <- cor(filtered_ME, filtered_cell_fraction, use='p', method = "spearman")
module_trait_cor_pvalue <- corPvalueStudent(module_trait_cor, nrow(filtered_ME))
module_trait_cor_adj_pvalue <- p.adjust(module_trait_cor_pvalue, method = 'BH')
module_trait_cor_adj_pvalue <- matrix(module_trait_cor_adj_pvalue, nrow = nrow(module_trait_cor_pvalue), ncol = ncol(module_trait_cor_pvalue))
colnames(module_trait_cor_adj_pvalue) <- colnames(module_trait_cor_pvalue)
rownames(module_trait_cor_adj_pvalue) <- rownames(module_trait_cor_pvalue)
write.table(module_trait_cor, 
            file = file.path(folder, "module_scRNAseq_cell_fraction_correlation_230720.txt"), 
            sep = '\t', quote = FALSE)
write.table(module_trait_cor_adj_pvalue, 
            file = file.path(folder, "module_scRNAseq_cell_fraction_correlation_P_230720.txt"), 
            sep = '\t', quote = FALSE)

#Edit rownames and colnames for output
rownames(module_trait_cor) <- gsub('ME', '', rownames(module_trait_cor))
rownames(module_trait_cor_adj_pvalue) <- gsub('ME', '', rownames(module_trait_cor_adj_pvalue))
module_trait_cor <- module_trait_cor[rownames(module_color_names), ]
module_trait_cor_adj_pvalue <- module_trait_cor_adj_pvalue[rownames(module_color_names), ]
rownames(module_trait_cor) <- as.character(module_color_names$Module_number)
rownames(module_trait_cor_adj_pvalue) <- as.character(module_color_names$Module_number)
module_trait_cor <- module_trait_cor[, c("Monocytes", "Granulocytes", "Erythroid.cells", "NK.cells", "Proliferative.Lymphocytes", "CD4..T.cells", "CD8..T.cells", "Dendritic.cells", 
                                         "B.cells", "Plasmablasts", "Cell.precursors..HSCs..MEPs..CMPs..GMPs.", "Megakaryocytes")]
module_trait_cor_adj_pvalue <- module_trait_cor_adj_pvalue[, c("Monocytes", "Granulocytes", "Erythroid.cells", "NK.cells", "Proliferative.Lymphocytes", "CD4..T.cells", "CD8..T.cells", "Dendritic.cells", 
                                                               "B.cells", "Plasmablasts", "Cell.precursors..HSCs..MEPs..CMPs..GMPs.", "Megakaryocytes")]
colnames(module_trait_cor) <- c("Monocytes", "Granulocytes", "Erythroid cells", "NK cells", "ProLymph", "CD4+ T cells", "CD8+ T cells", "DC", "B cells", "PB", "Precursors", "MK")

#Plot correlation heatmap as in figure 4B
pdf(file = file.path(folder, "module_scRNAseq_cell_fraction_correlation_250820.pdf"), 
    width=9, height = 8)
corrplot(module_trait_cor, method="color", type="full", col=viridis::viridis(30), p.mat = module_trait_cor_adj_pvalue, sig.level = c(.001, .01, .05), pch.cex=1.5, pch.col="white", 
         tl.col = "black", insig="label_sig", tl.srt = 45)
dev.off()

#Code for plotting module hub genes on single cell data as in figure S4-------
#Get sample list for scRNAseq samples
sample_list <- covid_info_file
sample_list <- subset(sample_list, sample_list$scRNAseqID != '')
sample_list$sample_name <- paste(sample_list$Name, sample_list$Pseudotime, sep = '_')

#Identify top 10 hub genes for each module
hub_genes <- data.frame()
for (i in 2:ncol(gene_module_membership)) {
  module_hub_genes <- as.character(gene_module_membership[order(gene_module_membership[,i], decreasing = TRUE), "Gene"][1:10])
  hub_genes <- rbind(hub_genes, data.frame(Gene=module_hub_genes, ModuleColor=gsub('ME', '', colnames(gene_module_membership)[i])))
}
hub_genes <- subset(hub_genes, hub_genes$ModuleColor != "grey")
hub_genes$Module <- module_color_to_number[as.character(hub_genes$ModuleColor), "Module_number"]
gene_list <- hub_genes$Gene

#Read average cell type data and filter for gene_list
cell_types <- c("Monocytes", "Granulocytes", "Erythroid", "NK", "Proliferative", "CD4", "CD8", "DC", "Bcells", "Plasmablasts", "Precursors", 
                "Megakaryocytes")
sample_names <- read.csv("output/Intergration_with_sc_data/Merge_sample_average.csv", row.names = 2)
folder <- "output/Intergration_with_sc_data/"
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

#Avergae expression for each pseudotime and cell type
pseudotime <- str_split_fixed(colnames(all_cell_type_data), '_', 5)[,4]
cell <- str_split_fixed(colnames(all_cell_type_data), '_', 5)[,5]
matDataHeat <- aggregate(t(all_cell_type_data), list(pseudotime, cell), mean)
rownames(matDataHeat) <- paste(matDataHeat$Group.2, matDataHeat$Group.1, sep = '_')
matDataHeat$Group.1 <- NULL
matDataHeat$Group.2 <- NULL

matDataHeat <- t(scale(matDataHeat))

#Assign modules to heatmap genes
matDataHeat_gene_module <- hub_genes[match(rownames(matDataHeat), hub_genes$Gene), c("Gene", "Module", "ModuleColor")]
matDataHeat_gene_module$Module <- factor(matDataHeat_gene_module$Module, levels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9","M10"))
matDataHeat_gene_module <- matDataHeat_gene_module[order(matDataHeat_gene_module$Module), ]
matDataHeat <- matDataHeat[as.character(matDataHeat_gene_module$Gene), ]

cell_colors <- c( "CD4" = "#EA811F", "CD8" = "maroon","NK" = "#6A3D8A",
                  "Bcells" = "goldenrod1", 'Plasmablasts' = 'goldenrod4', 
                  "Megakaryocytes" = "darkolivegreen",
                  'Precursors'='thistle1',
                  "DC" = "orchid3","Monocytes" = "paleturquoise3",  "Granulocytes" = "#1A3865",
                  "Erythroid" = "brown4", "Proliferative" = "#93D0AC")

col = as.character(unique(matDataHeat_gene_module$ModuleColor))
names(col) <- unique(matDataHeat_gene_module$Module)

#Plot heatmap as in figure S4
pseudotime <- str_split_fixed(colnames(matDataHeat), '_', 2)[,2]
cell <- str_split_fixed(colnames(matDataHeat), '_', 2)[,1]

col_ha = HeatmapAnnotation(df = data.frame(Celltype=cell, Pseudotime=pseudotime),
                           col = list(Celltype=cell_colors, 
                                      Pseudotime=c('0' = '#A7A9AC', '1' = '#E84F8C', '2' = '#7031B9', '3' = '#A04E9E', 
                                                   '4' = '#E65826', '5' = '#F99B1C', '6' = '#FDC077', '7' = '#51BBFE')), 
                           gp = gpar(col = "black"))
row_ha <- rowAnnotation(df = data.frame(Module = matDataHeat_gene_module$Module),
                        col = list(Module = col), 
                        gp = gpar(col = "black"))

pdf(file.path(folder, "Hub_genes_cell_type_expression.pdf"), width=10, height=8)
Heatmap(matDataHeat, col=colorRamp2(-5:5, rev(brewer.pal(11, "RdYlBu"))), 
        row_names_gp = gpar(fontsize = 6), left_annotation = row_ha, top_annotation = col_ha, 
        cluster_columns = FALSE, show_column_names = FALSE, 
        column_split = factor(cell, levels = cell_types), 
        column_title_gp = gpar(fontsize=0),
        row_split = matDataHeat_gene_module$Module, border = TRUE, cluster_row_slices = FALSE, row_title_rot = 0)
dev.off()


#Code to plot GSEA dot plot as in figure 4D---------
#Function to filter GSEA data to get upto n_terms significant terms
setwd(file.path(folder, "GSEA"))
get_module_gsea_data <- function(modules, n_terms){
  plist <- vector('list', length(modules))
  gsea_list <- c()
  
  for (i in 1:length(modules)) {
    module <- modules[i]
    file = list.files(path = file.path(module), pattern = "gsea_report_for_na_pos_\\d+\\.xls")
    gsea_results <- read.csv(file.path(module, file), sep = '\t')
    gsea_results <- subset(gsea_results, gsea_results$FDR.q.val < 0.05)
    gsea_results <- gsea_results[order(gsea_results$FDR.q.val), ]
    if (nrow(gsea_results) > n_terms) {
      gsea_results <- gsea_results[1:n_terms, ]
    }
    
    gsea_list <- c(gsea_list, as.character(gsea_results$NAME))
    
    plist[[i]] <- gsea_results
    
  }
  
  #Add p-values and NES value for each module
  p_vector <- c()
  nes_vector <- c()
  for (i in 1:length(modules)) {
    p_vector <- c(p_vector, paste("p", i, sep = '_'))
    nes_vector <- c(nes_vector, paste("nes", i, sep = '_'))
  }
  
  gsea_data <- data.frame(Term=gsea_list)
  gsea_data[, p_vector] <- NA
  gsea_data[, nes_vector] <- NA
  
  for (i in 1:length(modules)) {
    df <- plist[[i]]
    rownames(df) <- as.character(df$NAME)
    module_gsea_data <- df[as.character(gsea_data$Term),]
    gsea_data[, (i+1)] <- module_gsea_data$FDR.q.val
    gsea_data[, (i+1+length(modules))] <- module_gsea_data$NES
  }
  
  return(gsea_data)
}

#Select modules to plot GSEA results for
module_colors <- rev(c("bisque4", "tan", "brown", "orangered4", "darkturquoise", 
                       "plum1", "darkmagenta", "floralwhite", "lightpink3"))

#Get top 20 significant terms for each module
module_gsea_data <- get_module_gsea_data(module_colors, 20)
module_gsea_data <- module_gsea_data[order(module_gsea_data$nes_1, decreasing = FALSE), ]

#We selected a limited number of terms to plot
keep_terms <- read.table("keep_terms.txt")
keep_terms <- keep_terms$V1
module_selected_gsea_data <- subset(module_gsea_data, (module_gsea_data$Term %in% keep_terms)) 

#Reformat data 
p_vector <- colnames(module_gsea_data)[grep("p_", colnames(module_gsea_data))]
nes_vector <- colnames(module_gsea_data)[grep("nes_", colnames(module_gsea_data))]
selected_plot_data <- melt(module_selected_gsea_data[, c("Term", p_vector)], id.vars = "Term")
selected_plot_data$nes <- unlist(module_selected_gsea_data[, nes_vector])
selected_plot_data$Module <- module_colors[as.numeric(sub("p_", "", selected_plot_data$variable))]
selected_plot_data$value <- as.numeric(as.character(selected_plot_data$value))
selected_plot_data$p <- -1*log10(selected_plot_data$value+10^-6)
selected_plot_data$Term <- gsub('_', ' ', selected_plot_data$Term)
selected_plot_data$Term <- str_to_sentence(selected_plot_data$Term)
selected_plot_data$Term <- gsub('Go ', 'GO: ', selected_plot_data$Term)
selected_plot_data$Term <- gsub('Hallmark ', 'Hallmark: ', selected_plot_data$Term)
selected_plot_data$Term <- gsub('Kegg ', 'KEGG: ', selected_plot_data$Term)
selected_plot_data$Term2 <- reorder(selected_plot_data$Term, rev(1:nrow(selected_plot_data)))

module_names <- module_color_names[as.character(module_colors), "Module_number"]

#Plot GSEA dot plot as in figure 4D
pdf("Module_selected_GSEA_plot.pdf", width = 8, height = 8)
p <- ggplot(selected_plot_data, mapping = aes(x=Module, y=Term2, size=nes, color=p))
p <- p + geom_point()
p <- p + scale_x_discrete(limits=rev(c("bisque4", "tan", "brown", "orangered4", "darkturquoise", 
                                       "plum1", "darkmagenta", "floralwhite", "lightpink3")), 
                          labels=module_names)
p <- p + xlab("Module") + ylab("Term")
p <- p + scale_colour_gradient(high="#990000", low="#FF9999")
p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                            axis.text.x = element_text(size=14, color = "black", angle=45, hjust=1), 
                            axis.title = element_text(size = 20))
p
dev.off()
#Code for TFBS anaylsis and plotting results as in figure S3G-------
#Load the LOLA database
dbPath = "~/Projects/regions/LOLACore/hg38/"
regionDB = loadRegionDB(dbPath)

setwd(folder)

module_colors <- unique(gene_module$ModuleColor)

#Get the chromosomal locations of the promoter regions of the genes
gene_module_location <- inner_join(gene_module, unique_mart, by=c("Gene"="Gene_id"))
gene_module_location$promoter_start <- ifelse(gene_module_location$Strand == '+', gene_module_location$Start - 1500, gene_module_location$End - 500)
gene_module_location$promoter_end <- ifelse(gene_module_location$Strand == '+', gene_module_location$Start + 500, gene_module_location$End + 1500)

#Create background set
all_promoter_gr <- gene_module_location[, c("Chr", "promoter_start", "promoter_end", "Strand", "Gene")]
colnames(all_promoter_gr) <- c("Chromosome", "Start", "End","Strand", "gene_id")
all_promoter_gr <- GRanges(all_promoter_gr)

#Perform TFBS anaylsis for each module
all_tfbs_results <- data.frame()
top_hits <- c()

for (i in 1:length(module_colors)) {
  module_genes <- subset(gene_module_location, gene_module_location$ModuleColor == module_colors[i])
  
  #Create test set
  module_promoter_gr <- module_genes[, c("Chr", "promoter_start", "promoter_end", "Strand", "Gene")]
  colnames(module_promoter_gr) <- c("Chromosome", "Start", "End","Strand", "gene_id")
  module_promoter_gr <- GRanges(module_promoter_gr)
  
  #Run LOLA
  results <- runLOLA(module_promoter_gr, all_promoter_gr, regionDB, cores=1)
  
  write.table(results, paste("TFBS/", module_colors[i],"_promoter_LOLA_results.txt", sep = ''), 
              sep = '\t', quote = FALSE, row.names = FALSE)
  
  #Filter results for TFBS
  tfbs_results <- subset(results, results$collection == "encode_tfbs")
  tfbs_results$antibody <- str_split_fixed(tfbs_results$antibody, '_', 2)[,1]
  tfbs_results <- tfbs_results[!(duplicated(tfbs_results$antibody)==TRUE),]
  tfbs_results$Module <- module_colors[i]
  
  #Get top 5 significant hits
  all_tfbs_results <- rbind(all_tfbs_results, tfbs_results)
  tfbs_results <- subset(tfbs_results, tfbs_results$qValue < 0.05)
  top_hits <- c(top_hits, tfbs_results[order(tfbs_results$qValue), ]$antibody[1:5])
  
}
#Get all significant results
sig_all_tfbs_results <- subset(all_tfbs_results, all_tfbs_results$qValue < 0.05)

#Reformat data 
plot_data <- acast(sig_all_tfbs_results[, c("antibody", "Module", "oddsRatio")], antibody ~ Module)
plot_data <- as.data.frame(plot_data)
plot_data$grey <- NULL

top_hits <- unique(top_hits)
top_hits <- top_hits[!is.na(top_hits)]
plot_data <- plot_data[top_hits, ]

colnames(plot_data) <- module_color_names[colnames(plot_data), "Module_number"]
plot_data <- plot_data[, c("M1", "M2", "M3", "M5", "M6", "M7", "M10")]
plot_data <- plot_data[order(plot_data$M10, decreasing = TRUE), ]
plot_data <- plot_data[order(plot_data$M7, decreasing = TRUE), ]
plot_data <- plot_data[order(plot_data$M6, decreasing = TRUE), ]
plot_data <- plot_data[order(plot_data$M5, decreasing = TRUE), ]
plot_data <- plot_data[order(plot_data$M3, decreasing = TRUE), ]
plot_data <- plot_data[order(plot_data$M2, decreasing = TRUE), ]
plot_data <- plot_data[order(plot_data$M1, decreasing = TRUE), ]

#Plot heatmap as in figure S3G
pdf("TFBS/module_top5_TFBS.pdf", width = 4, height = 5)
Heatmap(plot_data, col=brewer.pal(9, "Reds")[2:9], cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, na_col = "white",
        row_split = 1:nrow(plot_data), column_split = 1:ncol(plot_data), row_gap = unit(0, "mm"), column_gap = unit(1, "mm"))
dev.off()


