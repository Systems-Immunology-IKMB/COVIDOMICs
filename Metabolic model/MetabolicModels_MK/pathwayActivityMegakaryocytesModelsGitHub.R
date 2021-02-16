#############################################################################
###
###   Differential pathway activity from FBA in Megakaryocytes
###    by Lena Best Nov. 2020, on R version 3.6.3 (2020-02-29)
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
df_MKsMetadata <- read.table(file = paste0(str_dataDir, "/MegakaryocytesMetadata.tsv"), 
                             header = T, sep = "\t", stringsAsFactors = F)

## read in active reaction counts per metabolic pathway
df_MKsPathwayAbundances <- read.table(file = paste0(str_dataDir, "/MegakaryocytesPathwayAbundances.tsv"), 
                                      header = T, row.names = 1, sep = "\t", stringsAsFactors = F)

#######################################################


###########################################################################
###
###   Part 1: Find differentially active pathways via Kruskal Wallis test
###
###########################################################################

#####################################
### Megakaryocytes
## Kruskal Wallis test
# Independent variables = pseudo time
df_independentVariables = df_MKsMetadata[,c( "CellID", "PseudoTime")]
df_independentVariables$PseudoTime <- factor(df_independentVariables$PseudoTime)
rownames(df_independentVariables) <- df_independentVariables$CellID
# dependent variable - matrix of active reactions per pathway
mtx_dependentVariable <- df_MKsPathwayAbundances[,rownames(df_independentVariables)]
# include only pathways that occur in at least any cell
mtx_dependentVariable <- mtx_dependentVariable[rowSums(mtx_dependentVariable) > 0,]
# rows and cols of both tables MUST match!
all.equal(rownames(df_independentVariables),colnames(mtx_dependentVariable))
# TRUE
# Run the test
df_kruskalPBMCpathwaysMK <- df_kruskalTest(mtx_dependentVariable, df_independentVariables$PseudoTime)
# [1] "Starting rank sum tests of 79 observations within 3912 samples! Mon Aug  3 17:46:27 2020"
# [1] "Finished all tests! Preparing results. Mon Aug  3 17:46:28 2020"
df_kruskalPBMCpathwaysMK <- df_kruskalPBMCpathwaysMK[order(abs(df_kruskalPBMCpathwaysMK$test.statistic), decreasing = T),]
nrow(df_kruskalPBMCpathwaysMK[df_kruskalPBMCpathwaysMK$p.adj <= 0.05,])
# 47
df_kruskalPBMCpathwaysMK[df_kruskalPBMCpathwaysMK$p.adj <= 0.05,]

#############################################################


###########################################################################
###
###   Part 2: Creating violin plots and heatmaps of pathways
###
###########################################################################

#############################################################
# create a violin plot of selected pathways:
# plot single result pathways as violins
str_result = "Pyruvate.metabolism"

# setup a postscript vector graphic
setEPS()
postscript(file = "vioplotMegakaryocytesPyruvateMetabolism.eps", 
           family = "ArialMT", fonts = "ArialMT",
           title = "Lena Best, Aug 2020", 
           bg = "transparent", fg = "black",
           width = 12, height = 9, # graphics dimension in inches
           paper = "special", pagecentre = T,
           colormodel = "srgb", useKerning = T, 
           pointsize = 10)
# adjust plot margins
par(mar=c(5,4,4,2)+.1)
# create a violinplot of active reactions in the selected pathway by pseudo time
vioplot::vioplot(as.numeric( mtx_dependentVariable[str_result,] ) ~ 
                 df_independentVariables[, c( "PseudoTime" ) ],
                 cex.axis = .7, las = 2, 
                 xlab = "Pseudotime", 
                 ylab = "Number of active reactions", 
                 main = str_result,
                 col = c("#A7A9AC", "#E84F8C", "#7031B9", "#A04E9E", "#E65826", "#F99B1C", "#FDC077", "#51BBFE"),
                 ylim = c(0, max(mtx_dependentVariable[str_result,]) + 1))
# add the number of single-cell models per pseudo time point
text(x = c(1:8), 
     y = rep( x = max(mtx_dependentVariable[str_result,]) + (max(mtx_dependentVariable[str_result,]) / 25), times = 8),
     labels = paste0("n = ", table( df_independentVariables[, c( "PseudoTime" ) ] )),
     cex = .6, col = "#808080")
# add mean-values as red bars
points(x = 1:8, 
       y = lst_dbl_meanByGroup( lst_numeric = as.numeric(mtx_dependentVariable[str_result,]), 
                                fct_groups = df_independentVariables[, c( "PseudoTime" ) ] ),
       pch = "_", col = 2, cex = 3)
dev.off()


#############################################################
## Create a Heatmap 
# define a p-value cutoff as significance threshold
dbl_pCutoff = 0.001  # (10^-3)
# create a matrix of active reactions per pathway
mtx_tmp <- matrix(ncol = length(levels(droplevels(df_independentVariables$PseudoTime))), nrow = 0)
# mean values are calculated for each pseudo time point, for each one of the top pathways
for(int_i in 1:sum(df_kruskalPBMCpathwaysMK$p.adj <= dbl_pCutoff)) {
  mtx_tmp <- rbind(mtx_tmp,
                   lst_dbl_meanByGroup( lst_numeric = as.numeric( mtx_dependentVariable[rownames(df_kruskalPBMCpathwaysMK)[df_kruskalPBMCpathwaysMK$p.adj <= dbl_pCutoff],][int_i,] ),
                                        fct_groups = droplevels( df_independentVariables$PseudoTime ) ))
}
# add names of the pathways as matrix rownames
rownames(mtx_tmp) <- rownames(df_kruskalPBMCpathwaysMK)[df_kruskalPBMCpathwaysMK$p.adj <= dbl_pCutoff]

# setup a postscript vector graphic
setEPS()
postscript(file = "MegakaryocytesHeatmap.eps", 
           family = "ArialMT", fonts = "ArialMT",
           title = "Lena Best, Aug 2020", 
           bg = "transparent", fg = "black",
           width = 8, height = 8, # graphics dimension in inches
           paper = "special", pagecentre = T,
           colormodel = "srgb", useKerning = T,
           pointsize = 10)
# adjust plot margins
par(mar = c(5,4,4,2)+.1)
# Create the actual heatmap plot
Heatmap(matrix = t(scale(t(mtx_tmp))), # scale the data
        col = colorRampPalette(brewer.pal(9, "Purples"))(100),
        color_space = "sRGB",
        name = "Scaled pathway activity",
        row_names_gp = gpar(cex = .8), 
        column_names_gp = gpar(cex = .9),
        column_names_rot = 0,
        column_names_centered = T, 
        column_names_side = "top",
        na_col = "grey", 
        cluster_rows = T, 
        show_row_dend = F,
        cluster_columns = F, 
        show_heatmap_legend = T, 
        column_title = "Megakaryocytes")
dev.off()

###########################################