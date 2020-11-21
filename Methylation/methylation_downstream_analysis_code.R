library(ggplot2)
library(reshape2)
library(UpSetR)
library(LOLA)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)


#Get number of differentially methylated sites (DMPs) between COVID pseudotimes and controls defined by the automatically selected rank cutoff--------
#Define pseudotimes and automatically selected rank cutoffs for each pseudotime
pseudotime <- c('1', '2', '3', '4', '5', '6')
auto_selected_rank_cutoff <- c(184451, 173594, 186970, 120294, 200838, 197180)
dmp_number_data <- data.frame()
listInput <- list()
hypo_list <- list()

#Read each differential methylation file one by one and extract the number of hypo and hyper methylated sites 
for (i in 1:length(pseudotime)) {
  folder <- paste("Control_vs_COVID_analysis/reports_Control_vs_COVID", pseudotime[i], sep='')
  dmp_data <- read.csv(file.path(folder, "differential_methylation_data/diffMethTable_site_cmp1.csv"))
  dmp_data <- subset(dmp_data, dmp_data$combinedRank <= auto_selected_rank_cutoff[i])
  hyper_len <- nrow(subset(dmp_data, dmp_data$mean.diff > 0))
  hypo_len <- nrow(subset(dmp_data, dmp_data$mean.diff < 0))
  dmp_data <- dmp_data[order(dmp_data$combinedRank), ]
  #Select top 30,000 DMPs for overlap comparison between pseudotimes
  listInput[[i]] <- as.character(dmp_data$cgid[1:30000])
  dmp_number_data <- rbind(dmp_number_data, data.frame(pseudotime[i], hyper_len, hypo_len))
  
}

dmp_number_data <- melt(dmp_number_data)

#Plot bar plot as in figure 4F
pdf(file = "Control_vs_COVID_analysis/dmp_number_comparison.pdf", width = 6, height = 7)
p <- ggplot(data = dmp_number_data, mapping = aes(x=pseudotime.i., y=value, fill=variable)) + 
  geom_bar(stat = "identity", position = "dodge", color="black")
p <- p + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
p <- p + xlab("Pseudotime") + ylab("Number of sites") + scale_x_discrete(labels = c("0 vs 1", "0 vs 2", "0 vs 3", "0 vs 4", "0 vs 5", "0 vs 6"))
p <- p + theme_bw() + theme(axis.text=element_text(size=14, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=12), 
                            legend.title = element_text(size=0), plot.title = element_text(size=20, face="bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
p
dev.off()

#Define intersections for overlap comparison between pseudotimes
names(listInput) <- paste("Pseudotime", pseudotime, sep = '')
intersections <- list(list("Pseudotime1"), list("Pseudotime2"), list("Pseudotime3"), list("Pseudotime4"), list("Pseudotime5"), list("Pseudotime6"),
                      list("Pseudotime1", "Pseudotime2"), list("Pseudotime2", "Pseudotime3"), 
                      list("Pseudotime1", "Pseudotime2", "Pseudotime3"), list("Pseudotime4", "Pseudotime5", "Pseudotime6"), 
                      list("Pseudotime3","Pseudotime4", "Pseudotime5", "Pseudotime6"), list("Pseudotime2", "Pseudotime3","Pseudotime4", "Pseudotime5", "Pseudotime6"),
                      list("Pseudotime1", "Pseudotime2", "Pseudotime3","Pseudotime4", "Pseudotime5", "Pseudotime6"))
#Plot intersection plot as in figure 4G
pdf("Control_vs_COVID_analysis/selected_intersections.pdf", width=8, height = 7)
upset(fromList(listInput), sets=names(listInput), nsets=6, keep.order = TRUE, nintersects = NA, 
      mainbar.y.label = "DMP Intersections", sets.x.label = "#DMPs", intersections = intersections,
      text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3), number.angles = 30)
dev.off()

#Get number of differentially methylated sites (DMPs) between subsequent COVID pseudotimes defined by the automatically selected rank cutoff--------
#Define comparisons and automatically selected rank cutoffs for each pseudotime

comparisons <- c("1_vs_3", "2_vs_3","3_vs_4", "4_vs_5", "5_vs_6")

auto_selected_rank_cutoff <- c(0,0, 80974, 137377, 196182)
dmp_number_data <- data.frame()
listInput <- list()

#Read each differential methylation file one by one and extract the number of hypo and hyper methylated sites 
for (i in 1:length(comparisons)) {
  folder <- paste("COVID_longitudinal_analysis/reports", comparisons[i], sep='_')
  dmp_data <- read.csv(file.path(folder, "differential_methylation_data/diffMethTable_site_cmp1.csv"))
  dmp_data <- subset(dmp_data, dmp_data$combinedRank <= auto_selected_rank_cutoff[i])
  hyper_len <- nrow(subset(dmp_data, dmp_data$mean.diff < 0))
  hypo_len <- nrow(subset(dmp_data, dmp_data$mean.diff > 0))
  dmp_data <- dmp_data[order(dmp_data$combinedRank), ]
  listInput[[i]] <- as.character(dmp_data$cgid)
  dmp_number_data <- rbind(dmp_number_data, data.frame(comparisons[i], hyper_len, hypo_len))
}

dmp_number_data <- melt(dmp_number_data)

#Plot bar plot as in figure S5B
pdf(file = "COVID_longitudinal_analysis/dmp_number_comparison.pdf", width = 5.5, height = 6)
p <- ggplot(data = dmp_number_data, mapping = aes(x=comparisons.i., y=value, fill=variable)) + 
  geom_bar(stat = "identity", position = "dodge", color="black")
p <- p + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
p <- p + xlab("Comparisons") + ylab("Number of sites") + scale_x_discrete(labels = c("1 vs 3", "2 vs 3", "3 vs 4", "4 vs 5", "5 vs 6"))
p <- p + theme_bw() + theme(axis.text=element_text(size=14, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=12), 
                            legend.title = element_text(size=0), plot.title = element_text(size=20, face="bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
p
dev.off()

#LOLA analysis for TFBS enrichment analysis-----------

#Load region data base
dbPath = "../../../regions/LOLACore/hg19/"
regionDB = loadRegionDB(dbPath)

#Define pseudotimes and automatically selected rank cutoffs for each pseudotime
pseudotime <- c('1', '2', '3', '4', '5', '6')
auto_selected_rank_cutoff <- c(184451, 173594, 186970, 120294, 200838, 197180)

for (i in 1:length(pseudotime)) {
  #Read DMP file
  folder <- paste("Control_vs_COVID_analysis/reports_Control_vs_COVID", pseudotime[i], sep='')
  all_sites <- read.csv(file.path(folder, "differential_methylation_data/diffMethTable_site_cmp1.csv"))
  
  #Extract chromosomal locations of hypo and hypermethylated sites
  dmp_sites <- subset(all_sites, all_sites$combinedRank <= auto_selected_rank_cutoff[i])
  hypo_sites <- subset(dmp_sites, dmp_sites$mean.diff < 0)
  hyper_sites <- subset(dmp_sites, dmp_sites$mean.diff > 0)
  
  hypo_data <- hypo_sites[, c("Chromosome", "Start", "Strand", "cgid")]
  hypo_data$End <- hypo_data$Start + 1
  hypo_gr <- GRanges(seqnames = hypo_data$Chromosome,
                     ranges = IRanges(hypo_data$Start, hypo_data$End),
                     strand = hypo_data$Strand, 
                     id = hypo_data$cgid)
  
  hyper_data <- hyper_sites[, c("Chromosome", "Start", "Strand")]
  hyper_data$End <- hyper_data$Start + 1
  hyper_gr <- GRanges(seqnames = hyper_data$Chromosome,
                      ranges = IRanges(hyper_data$Start, hyper_data$End),
                      strand = hyper_data$Strand)
  
  #Set background set
  universe_data <- all_sites[, c("Chromosome", "Start", "Strand")]
  universe_data$End <- universe_data$Start + 1
  universe_gr <- GRanges(seqnames = universe_data$Chromosome,
                         ranges = IRanges(universe_data$Start, universe_data$End),
                         strand = universe_data$Strand)
  
  #Run LOLA and output results
  hypo_results <- runLOLA(hypo_gr, universe_gr, regionDB, cores=1)
  write.table(hypo_results, file.path(folder, paste("Pseudotime", pseudotime[i],"DMP_hypo_LOLA_results.txt", sep = '_')), 
              sep = '\t', quote = FALSE, row.names = FALSE)
  
  hyper_results <- runLOLA(hyper_gr, universe_gr, regionDB, cores=1)
  write.table(hyper_results, file.path(folder, paste("Pseudotime", pseudotime[i],"DMP_hyper_LOLA_results.txt", sep = '_')), 
              sep = '\t', quote = FALSE, row.names = FALSE)
  
  #Extract results for transcription factor binding sites (TFBS)
  hypo_tfbs_results <- subset(hypo_results, hypo_results$collection == "encode_tfbs")
  
  hyper_tfbs_results <- subset(hyper_results, hyper_results$collection == "encode_tfbs")
  
  
  hypo_tfbs_results$direction <- "Hypomethylated"
  hyper_tfbs_results$direction <- "Hypermethylated"
  
  hypo_tfbs_results$antibody <- str_split_fixed(hypo_tfbs_results$antibody, '_', 2)[,1]
  hyper_tfbs_results$antibody <- str_split_fixed(hyper_tfbs_results$antibody, '_', 2)[,1]
  
  hypo_tfbs_results <- hypo_tfbs_results[!(duplicated(hypo_tfbs_results$antibody)==TRUE),]
  hyper_tfbs_results <- hyper_tfbs_results[!(duplicated(hyper_tfbs_results$antibody)==TRUE),]
  
  #Filter for TFBS significantly enriched in hyper and/or hypo methylated sites
  tfbs_results <- full_join(hypo_tfbs_results, hyper_tfbs_results, by="antibody", suffix=c(".hypo", ".hyper"))
  tfbs_results <- subset(tfbs_results, tfbs_results$qValue.hypo < 0.05 | tfbs_results$qValue.hyper < 0.05)
  
  #Select top 5 TFBS for plotting
  top_hypo <- tfbs_results[order(tfbs_results$qValue.hypo), ]$antibody[1:5]
  top_hyper <- tfbs_results[order(tfbs_results$qValue.hyper), ]$antibody[1:5]
  
  #Plot results as heatmap as in figure 4H
  plot_data <- subset(tfbs_results, tfbs_results$antibody %in% c(top_hypo, top_hyper))
  plot_data <- plot_data[, c("oddsRatio.hypo", "oddsRatio.hyper", "antibody")]
  plot_data <- as.data.frame(plot_data)
  rownames(plot_data) <- as.character(plot_data$antibody)
  plot_data$antibody <- NULL
  colnames(plot_data) <- c("Hypo", "Hyper")
  
  pdf(file.path(folder, paste("Pseudotime", pseudotime[i],"DMP_LOLA_heatmap.pdf", sep = '_')), width = 3, height = 4)
  print(Heatmap(plot_data, col=brewer.pal(9, "Reds"), cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, 
                row_split = 1:nrow(plot_data), column_split = 1:ncol(plot_data), row_gap = unit(0, "mm"), column_gap = unit(0, "mm")))
  dev.off()
}

#Code for cellular deconvolution using methylation data-------
rnb.path <- "Analysis/reports_preprocessing/rnbSet_preprocessed/"#Specify the path to your preprocessed RnBSet
ref.path <- "Analysis/Reinius_Blood_Reference/"#Path to the Reinius reference
plot.path <- "Analysis/reports_cell_estimation/cell_estimation_figures/"#Some path to store a plot overviewing the results
new.path <- "Analysis/reports_cell_estimation/cell_estimation_rnbset/"#Path to store the RnBSet with contributions added to the phenotypic information

#Load RnBeads objects for test and reference sets
rnb.set <- load.rnb.set(rnb.path)
ref.set <- load.rnb.set(ref.path)

#Edit phenotype data for reference set
ph <- pheno(ref.set)
rem.samples <- ph$'Cell/Tissue' %in% c("PBMC","granulocytes","whole blood")
ref.set <- remove.samples(ref.set,rem.samples)

#Edit phenotype data for test set
rnb.set2 <- addPheno(rnb.set, pheno(rnb.set)[, "Sentrix_ID"], "Sentrix ID")
rnb.set2 <- addPheno(rnb.set2, pheno(rnb.set)[, "Sentrix_Position"], "Sentrix Position")
rnb.set2 <- addPheno(rnb.set2, pheno(rnb.set)[, "Sample_ID"], "ID")
rnb.set3 <- rnb.combine.arrays(rnb.set2,ref.set,type="all.x")
rnb.set <- rnb.set3
rnb.set <- rnb.execute.imputation(rnb.set)

#Run RnBeads for cell type estimation
rnb.options(
  "inference.reference.methylome.column"="Cell Type",
  #' number of most variable sites to be selected for potential link to the reference cell types
  "inference.max.cell.type.markers"=50000,
  #' number of sites used for the final calculation of cell type contributions
  "inference.top.cell.type.markers"=500
)
cell.est <- rnb.execute.ct.estimation(rnb.set,cell.type.column="Cell/Tissue")

#Add cell type estimations to original RnBeads object and save
add.pheno <- cell.est$contributions.nonneg
rnb.set <- load.rnb.set(rnb.path)
for(i in 1:ncol(add.pheno)){
  rnb.set <- addPheno(rnb.set,add.pheno[,i],colnames(add.pheno)[i])
}
save.rnb.set(rnb.set,file.path(new.path,"set_with_contributions"))
write.table(pheno(rnb.set), file.path(plot.path, "contributions.txt"), quote = FALSE, sep = '\t')
epi_cell_fraction <- pheno(rnb.set)

#Combine column sets
epi_cell_fraction$Granulocytes <- epi_cell_fraction$eosinophils + epi_cell_fraction$neutrophils

#Rename
colnames(epi_cell_fraction)[which(names(epi_cell_fraction) == "CD56..NK.cells")] <- "NK cells"
colnames(epi_cell_fraction)[which(names(epi_cell_fraction) == "CD19..B.cells")] <- "B cells"
colnames(epi_cell_fraction)[which(names(epi_cell_fraction) == "CD14..monocytes")] <- "Monocytes"
colnames(epi_cell_fraction)[which(names(epi_cell_fraction) == "CD4..T.cells")] <- "CD4+ T cells"
colnames(epi_cell_fraction)[which(names(epi_cell_fraction) == "CD8..T.cells")] <- "CD8+ T cells"
colnames(epi_cell_fraction)[which(names(epi_cell_fraction) == "Sample_Name2")] <- "Sample"

#subset column
epi_cell_fraction <- subset(epi_cell_fraction, select = c("Sample_Name", "B cells","CD4+ T cells", "CD8+ T cells", "Granulocytes", "Monocytes", 
                                                          "NK cells"))

# pivot longer
epi_cell_fraction_pivot <- pivot_longer(epi_cell_fraction, cols = 2:7, names_to = "Cell_Type")
colnames(epi_cell_fraction_pivot)[which(names(epi_cell_fraction_pivot) == "value")] <- "Contribution" 

#Plot stacked bar plot as in figure S5A
pdf("RnBeads_contribution_controls_incl.pdf", width=6, height = 6)
p <- ggplot(epi_cell_fraction_pivot, aes(x = Sample, y = Contribution, fill = Cell_Type)) + 
  geom_bar(position="fill", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()

p <- p + scale_fill_manual(values=c( "CD4+ T cells" = "#EA811F", "CD8+ T cells" = "maroon","NK cells" = "#6A3D8A",
                                     "B cells" = "goldenrod1", "Monocytes" = "paleturquoise3",  "Granulocytes" = "#1A3865"))

p <- p + theme_bw() + theme(axis.text = element_text(size = 14, colour = "black"), axis.title = element_text(size = 16))
p
dev.off()



