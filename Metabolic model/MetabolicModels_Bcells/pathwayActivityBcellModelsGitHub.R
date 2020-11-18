#############################################################################
###
###   Differential pathway activity from Flux Balance Analysis in B-cells
###          by Lena Best Nov. 2020, on R version 3.6.3 (2020-02-29)
###
#############################################################################

#############################################################
library(coin)
library(parallel)
library(RColorBrewer)
library(ComplexHeatmap)

# path to where the data-files are stored (defaults to current working directory)
str_dataDir <- getwd()


#############################################################
## define a function to run kruskal-wallis test in parallel on a data.frame or matrix of reaction counts
df_kruskalTest <- function(df_dependantVariable, lst_fct_metaData, int_cores = 7) {
  require("coin")
  require(parallel)
  
  ## make sure metaData is a factor
  lst_fct_metaData <- droplevels(as.factor(lst_fct_metaData))

  ## make sure dependant data is numeric
  mtx_dependantVariable <- data.matrix( frame = df_dependantVariable )
  
  ## retain rownames
  # add a column for the rownames (can later be used to get them back)
  # this column needs to be numeric as well in order to be succesfully passed to a function!
  # thus we employ factors
  int_colNum <- ncol(mtx_dependantVariable)
  fct_rowNames <- factor(x = rownames(df_dependantVariable))
  mtx_dependantVariable <- cbind(mtx_dependantVariable, as.numeric( fct_rowNames ))
  
  #Internal function for calling multiple correlations in parrallel mode
  df_internalRankTest <- function(row_dependantVariable) {
    int_rowName <- row_dependantVariable[int_colNum + 1]
    row_dependantVariable <- as.numeric( row_dependantVariable[1 : int_colNum] )
    
    flt_pValue <- 2
    flt_HValue <- 0
    
    # run the kruskal test
    #H0: row_dependantVariable (e.g. Expression of Gene Xyz) is independent of lst_var_metaData (e.g. Age)
    #H1: row_dependantVariable (e.g. Expression of Gene Xyz) is altered by the samples lst_var_metaData (e.g. Age)
    try( 
      obj_kruskalTest <- kruskal_test( row_dependantVariable ~ lst_fct_metaData, 
                                      alternative = "two.sided", 
                                      paired = F, 
                                      exact = T, 
                                      ties.method = "average-scores", 
                                      correct = F)
    )
    if( exists("obj_kruskalTest") ) {
      flt_HValue <- obj_kruskalTest@statistic@teststatistic
      flt_pValue <- obj_kruskalTest@distribution@pvalue( flt_HValue )
    }
    else print(paste("Unable to test dependant variable row",fct_rowNames[int_rowName],"!"))
    
    return( c(int_rowName, flt_pValue, flt_HValue) )
  }
  
  # start working in parrallel
  obj_clusters <- makeCluster(int_cores, type="SOCK")
  clusterCall(cl = obj_clusters, fun = function() library(coin))
  clusterExport(cl = obj_clusters, varlist = c("mtx_dependantVariable","int_colNum","lst_fct_metaData","fct_rowNames","df_internalRankTest"), envir = environment())
  
  #info to user
  print( paste( "Starting Kruskal Wallis tests of",nrow(df_dependantVariable),"observations within",int_colNum,"samples!", date() ))
  
  # run rank sum tests in parrallel
  df_resultsKruskal <- as.data.frame( t( parApply(cl = obj_clusters, X = mtx_dependantVariable, MARGIN = 1, FUN = df_internalRankTest) ))
  # df_resultsKruskal <- as.data.frame( t( apply(X = mtx_dependantVariable, MARGIN = 1, FUN = df_internalRankTest) ))

  # clean up
  stopCluster(obj_clusters)
  
  #info to user
  print( paste("Finished all tests! Preparing results.", date()) )
  
  colnames(df_resultsKruskal) <- c("n","p.value","test.statistic")
  
  # correct p-values for multiple testing
  df_resultsKruskal$p.adj <- p.adjust(p = df_resultsKruskal$p.value, method = "BH")
  
  return( df_resultsKruskal )
}

## function to return mean values from a list of numeric values, grouped by a factor
lst_dbl_meanByGroup <- function(lst_numeric, fct_groups, na.rm = T) {
  fct_groups <- droplevels( as.factor(fct_groups) )

  lst_dbl_result <- NULL
  for(int_lvl in levels(fct_groups)) {
    lst_dbl_result <- c(lst_dbl_result, mean(lst_numeric[fct_groups == int_lvl], na.rm = na.rm) )
  }
  names(lst_dbl_result) <- levels(fct_groups)
  return (lst_dbl_result)
}


#############################################################
## read in metadata
df_BcellsMetadata <- read.table(file = paste0(str_dataDir, "/BcellsMetadata.tsv"), 
                                header = T, sep = "\t", stringsAsFactors = F)

## read in active reaction counts per metabolic pathway
df_BcellsPathwayAbundances <- read.table(file = paste0(str_dataDir, "/BcellsPathwayAbundances.tsv"), 
                                         header = T, row.names = 1, sep = "\t", stringsAsFactors = F)

#######################################################


###########################################################################
###
###   Part 1: Find differentially active pathways via Kruskal Wallis test
###
###########################################################################

#####################################
### Plasmablasts
## Kruskal Wallis test
# Independent variables = pseudo time
df_independentVariables = df_BcellsMetadata[df_BcellsMetadata$CellType == "Plasmablast",
                                            c( "CellID", "PseudoTime")]
df_independentVariables$PseudoTime <- factor(df_independentVariables$PseudoTime)
rownames(df_independentVariables) <- df_independentVariables$CellID
# dependent variable - matrix of active reactions per pathway
mtx_dependentVariable <- df_BcellsPathwayAbundances[,rownames(df_independentVariables)]
# include only pathways that occur in at least any cell
mtx_dependentVariable <- mtx_dependentVariable[rowSums(mtx_dependentVariable) > 0,]
# rows and cols of both tables MUST match!
all.equal(rownames(df_independentVariables),colnames(mtx_dependentVariable))
# TRUE
# Run the test
df_kruskalPBMCpathwaysPlasmablasts <- df_kruskalTest(mtx_dependentVariable, df_independentVariables$PseudoTime)
# [1] "Starting Kruskal Wallis tests of 81 observations within 5539 samples! Wed Nov 18 15:27:02 2020"
# [1] "Finished all tests! Preparing results. Wed Nov 18 15:27:03 2020"
df_kruskalPBMCpathwaysPlasmablasts <- df_kruskalPBMCpathwaysPlasmablasts[order(abs(df_kruskalPBMCpathwaysPlasmablasts$test.statistic), decreasing = T),]
nrow(df_kruskalPBMCpathwaysPlasmablasts[df_kruskalPBMCpathwaysPlasmablasts$p.adj <= 0.05,])
# 71
df_kruskalPBMCpathwaysPlasmablasts[df_kruskalPBMCpathwaysPlasmablasts$p.adj <= 0.05,]


#############################################################
### B-cells (Memory subgroup)
## Kruskal Wallis test
# Independent variables = pseudo time
df_independentVariables = df_BcellsMetadata[df_BcellsMetadata$CellType == "Memory_B-cell",
                                            c( "CellID", "PseudoTime")]
df_independentVariables$PseudoTime <- factor(df_independentVariables$PseudoTime)
rownames(df_independentVariables) <- df_independentVariables$CellID
# dependent variable - matrix of active reactions per pathway
mtx_dependentVariable <- df_BcellsPathwayAbundances[,rownames(df_independentVariables)]
# include only pathways that occur in at least any cell
mtx_dependentVariable <- mtx_dependentVariable[rowSums(mtx_dependentVariable) > 0,]
# rows and cols of both tables MUST match!
all.equal(rownames(df_independentVariables),colnames(mtx_dependentVariable))
# TRUE
# Run the test
df_kruskalPBMCpathwaysMemoryBcells <- df_kruskalTest(mtx_dependentVariable, df_independentVariables$PseudoTime )
# [1] "Starting Kruskal Wallis tests of 81 observations within 10894 samples! Wed Nov 18 15:35:04 2020"
# [1] "Finished all tests! Preparing results. Wed Nov 18 15:35:05 2020"
df_kruskalPBMCpathwaysMemoryBcells <- df_kruskalPBMCpathwaysMemoryBcells[order(abs(df_kruskalPBMCpathwaysMemoryBcells$test.statistic), decreasing = T),]
df_kruskalPBMCpathwaysMemoryBcells[df_kruskalPBMCpathwaysMemoryBcells$p.adj <= 0.05,]
nrow(df_kruskalPBMCpathwaysMemoryBcells[df_kruskalPBMCpathwaysMemoryBcells$p.adj <= 0.05,])
# 78


#############################################################
### B-cells (Naiive subgroup)
## Kruskal Wallis test
# Independent variables = pseudo time
df_independentVariables = df_BcellsMetadata[df_BcellsMetadata$CellType == "Naiive_B-cell",
                                            c( "CellID", "PseudoTime")]
df_independentVariables$PseudoTime <- factor(df_independentVariables$PseudoTime)
rownames(df_independentVariables) <- df_independentVariables$CellID
# dependent variable - matrix of active reactions per pathway
mtx_dependentVariable <- df_BcellsPathwayAbundances[,rownames(df_independentVariables)]
# include only pathways that occur in at least any cell
mtx_dependentVariable <- mtx_dependentVariable[rowSums(mtx_dependentVariable) > 0,]
# rows and cols of both tables MUST match!
all.equal(rownames(df_independentVariables),colnames(mtx_dependentVariable))
# TRUE
# Run the test
df_kruskalPBMCpathwaysNaiiveBcells <- df_kruskalTest(mtx_dependentVariable, df_independentVariables$PseudoTime )
# [1] "Starting Kruskal Wallis tests of 80 observations within 3623 samples! Wed Nov 18 15:35:48 2020"
# [1] "Finished all tests! Preparing results. Wed Nov 18 15:35:48 2020"
df_kruskalPBMCpathwaysNaiiveBcells <- df_kruskalPBMCpathwaysNaiiveBcells[order(abs(df_kruskalPBMCpathwaysNaiiveBcells$test.statistic), decreasing = T),]
df_kruskalPBMCpathwaysNaiiveBcells[df_kruskalPBMCpathwaysNaiiveBcells$p.adj <= 0.05,]
nrow(df_kruskalPBMCpathwaysNaiiveBcells[df_kruskalPBMCpathwaysNaiiveBcells$p.adj <= 0.05,])
# 73

#############################################################


###########################################################################
###
###   Part 2: Select a list of candidate pathways differential
###             abundant in all 3 cell sub types
###
###########################################################################

#############################################################
## Create a merged table of p-values from all three cell types
df_kruskalPBMCpathwaysAllBcells <- merge(x = df_kruskalPBMCpathwaysPlasmablasts[,"p.adj",drop = F], 
                                         y = df_kruskalPBMCpathwaysMemoryBcells[,"p.adj",drop = F], 
                                         by = "row.names", all = T)
rownames(df_kruskalPBMCpathwaysAllBcells) <- df_kruskalPBMCpathwaysAllBcells$Row.names
df_kruskalPBMCpathwaysAllBcells$Row.names <- NULL
colnames(df_kruskalPBMCpathwaysAllBcells) <- c("Plasmablasts", "MemoryBcells")
#
df_kruskalPBMCpathwaysAllBcells <- merge(x = df_kruskalPBMCpathwaysAllBcells,
                                         y = df_kruskalPBMCpathwaysNaiiveBcells[,"p.adj",drop = F], 
                                         by = "row.names", all = T)
rownames(df_kruskalPBMCpathwaysAllBcells) <- df_kruskalPBMCpathwaysAllBcells$Row.names
df_kruskalPBMCpathwaysAllBcells$Row.names <- NULL
colnames(df_kruskalPBMCpathwaysAllBcells) <- c("Plasmablasts", "MemoryBcells","NaiiveBcells")

# mask missing values (which occur from merging tables) with 2 (real p-values will only range from 0 to 1)
df_kruskalPBMCpathwaysAllBcells[is.na(df_kruskalPBMCpathwaysAllBcells)] <- 2

# order by p-values of plasmablasts (or alternatively by memory b-cells and last by naiive b-cells)
df_kruskalPBMCpathwaysAllBcells <- df_kruskalPBMCpathwaysAllBcells[order(df_kruskalPBMCpathwaysAllBcells$Plasmablasts, 
                                                                         df_kruskalPBMCpathwaysAllBcells$MemoryBcells, 
                                                                         df_kruskalPBMCpathwaysAllBcells$NaiiveBcells, 
                                                                         decreasing = F),]
# define a p-value cutoff as significance threshold
dbl_pCutoff = 0.001  # (10^-3)
# names of pathways beeing significant in all 3 datasets (... >= 3)
lst_str_topPathways <- rownames(df_kruskalPBMCpathwaysAllBcells)[( ((df_kruskalPBMCpathwaysAllBcells$Plasmablasts <= dbl_pCutoff) + 0) + 
                                                                   ((df_kruskalPBMCpathwaysAllBcells$MemoryBcells <= dbl_pCutoff) + 0) + 
                                                                   ((df_kruskalPBMCpathwaysAllBcells$NaiiveBcells <= dbl_pCutoff) + 0)) >= 3]
length(lst_str_topPathways)
# 41
# pathway called "Miscellaneous" will be dropped
lst_str_topPathways <- lst_str_topPathways[lst_str_topPathways != "Miscellaneous"]
# 40
# limit to only top 20
lst_str_topPathways <- lst_str_topPathways[1:20]

#############################################################


###########################################################################
###
###   Part 3: Based of the list of top pathways create a combined heatmap
###
###########################################################################

#####################################
### Plasmablasts
df_independentVariables = df_BcellsMetadata[df_BcellsMetadata$CellType == "Plasmablast",
                                            c( "CellID", "PseudoTime")]
df_independentVariables$PseudoTime <- factor(df_independentVariables$PseudoTime)
rownames(df_independentVariables) <- df_independentVariables$CellID
# dependent variable - matrix of active reactions per pathway
mtx_dependentVariable <- df_BcellsPathwayAbundances[,rownames(df_independentVariables)]
# include only pathways that occur in at least any cell
mtx_dependentVariable <- mtx_dependentVariable[rowSums(mtx_dependentVariable) > 0,]
# rows and cols of both tables MUST match!
all.equal(rownames(df_independentVariables),colnames(mtx_dependentVariable))
# TRUE
# create a matrix of active reactions per pathway
mtx_tmp <- matrix(ncol = length(levels(df_independentVariables$PseudoTime)), nrow = 0)
# mean values are calculated for each pseudo time point, for each one of the top pathways
for(int_i in 1:length(lst_str_topPathways)) {
  mtx_tmp <- rbind(mtx_tmp,
                   lst_dbl_meanByGroup(lst_numeric = as.numeric(mtx_dependentVariable[lst_str_topPathways[int_i],]),
                                       fct_groups = df_independentVariables$PseudoTime))
}
rownames(mtx_tmp) <- lst_str_topPathways
# rename the temporary matrix
mtx_Plasmablasts <- mtx_tmp; rm(mtx_tmp, df_independentVariables, mtx_dependentVariable)

#####################################
### B-cells (Memory subgroup)
df_independentVariables = df_BcellsMetadata[df_BcellsMetadata$CellType == "Memory_B-cell",
                                            c( "CellID", "PseudoTime")]
df_independentVariables$PseudoTime <- factor(df_independentVariables$PseudoTime)
rownames(df_independentVariables) <- df_independentVariables$CellID
# dependent variable - matrix of active reactions per pathway
mtx_dependentVariable <- df_BcellsPathwayAbundances[,rownames(df_independentVariables)]
# include only pathways that occur in at least any cell
mtx_dependentVariable <- mtx_dependentVariable[rowSums(mtx_dependentVariable) > 0,]
# rows and cols of both tables MUST match!
all.equal(rownames(df_independentVariables),colnames(mtx_dependentVariable))
# TRUE
# create a matrix of active reactions per pathway
mtx_tmp <- matrix(ncol = length(levels(df_independentVariables$PseudoTime)), nrow = 0)
# mean values are calculated for each pseudo time point, for each one of the top pathways
for(int_i in 1:length(lst_str_topPathways)) {
  mtx_tmp <- rbind(mtx_tmp,
                   lst_dbl_meanByGroup(lst_numeric = as.numeric(mtx_dependentVariable[lst_str_topPathways[int_i],]),
                                       fct_groups = df_independentVariables$PseudoTime))
}
rownames(mtx_tmp) <- lst_str_topPathways
# rename the temporary matrix
mtx_MemoryBcells <- mtx_tmp; rm(mtx_tmp, df_independentVariables, mtx_dependentVariable)

#####################################
### B-cells (Naiive subgroup)
df_independentVariables = df_BcellsMetadata[df_BcellsMetadata$CellType == "Naiive_B-cell",
                                            c( "CellID", "PseudoTime")]
df_independentVariables$PseudoTime <- factor(df_independentVariables$PseudoTime)
rownames(df_independentVariables) <- df_independentVariables$CellID
# dependent variable - matrix of active reactions per pathway
mtx_dependentVariable <- df_BcellsPathwayAbundances[,rownames(df_independentVariables)]
# include only pathways that occur in at least any cell
mtx_dependentVariable <- mtx_dependentVariable[rowSums(mtx_dependentVariable) > 0,]
# rows and cols of both tables MUST match!
all.equal(rownames(df_independentVariables),colnames(mtx_dependentVariable))
# TRUE
# create a matrix of active reactions per pathway
mtx_tmp <- matrix(ncol = length(levels(df_independentVariables$PseudoTime)), nrow = 0)
# mean values are calculated for each pseudo time point, for each one of the top pathways
for(int_i in 1:length(lst_str_topPathways)) {
  mtx_tmp <- rbind(mtx_tmp,
                   lst_dbl_meanByGroup(lst_numeric = as.numeric(mtx_dependentVariable[lst_str_topPathways[int_i],]),
                                       fct_groups = df_independentVariables$PseudoTime))
}
rownames(mtx_tmp) <- lst_str_topPathways
# rename the temporary matrix
mtx_NaiiveBcells <- mtx_tmp; rm(mtx_tmp, df_independentVariables, mtx_dependentVariable)


#####################################
## Create a Heatmap with all three cell types combined
setEPS()
postscript(file = "PlasmablastsBcellsCombinedHeatmap.eps", 
           family = "ArialMT", fonts = "ArialMT",
           title = "Lena Best, Aug 2020", 
           bg = "transparent", fg = "black",
           width = 10, height = 5, # graphics dimension in inches
           paper = "special", pagecentre = T,
           colormodel = "srgb", useKerning = T, #"cmyk"
           pointsize = 10)
Heatmap(matrix = t( scale( t(mtx_Plasmablasts))),
        col = colorRampPalette(brewer.pal(9, "Purples"))(100),
        color_space = "sRGB",
        name = "Plasmablasts", #"Scaled pathway activity",
        row_names_gp = gpar(cex = .8), 
        column_names_gp = gpar(cex = .9),
        column_names_rot = 0,
        column_names_centered = T, 
        column_names_side = "top",
        na_col = "grey", 
        cluster_rows = T, 
        show_row_dend = F,
        cluster_columns = F, 
        show_heatmap_legend = F, 
        column_title = "Plasmablasts") +
Heatmap(matrix = t (scale( t(mtx_MemoryBcells))), 
        col = colorRampPalette(brewer.pal(9, "Purples"))(100),
        color_space = "sRGB",
        name = "Scaled pathway activity", #"Memory B-cells",
        row_names_gp = gpar(cex = .8),
        column_names_gp = gpar(cex = .9),
        column_names_rot = 0,
        column_names_centered = T,
        column_names_side = "top",
        na_col = "grey", 
        cluster_rows = F, 
        show_row_dend = F,
        cluster_columns = F, 
        show_heatmap_legend = T,
        column_title = "Memory B-cells") + 
Heatmap(matrix = t( scale( t(mtx_NaiiveBcells))), 
        col = colorRampPalette(brewer.pal(9, "Purples"))(100),
        color_space = "sRGB",
        name = "Naiive B-cells",
        row_names_gp = gpar(cex = .8),
        column_names_gp = gpar(cex = .9),
        column_names_rot = 0, 
        column_names_centered = T,
        column_names_side = "top",
        na_col = "grey", 
        cluster_rows = F, 
        cluster_columns = F, 
        show_heatmap_legend = F,
        column_title = "Naiive B-cells")
dev.off()

###########################################