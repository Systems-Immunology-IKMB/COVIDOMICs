#
# @author:	Axel Künstner
# @email:	axel.kuenstner@uni-luebeck.de
# @date:	2020-08-13
#
# Dorothea/CARNIVAL workflow was done according to
# https://github.com/saezlab/NetworkModeling_course_2020
#

library(CARNIVAL)

library(tidyverse)
library(curl)
library(viper)
library(OmnipathR)


# Get data ----------------------------------------------------------------

cova            <- readRDS("../data/covariates.RDS")
log.norm.counts <- readRDS("../data/RNAseq_cnts_norm.RDS")

cova$Pseudotime <- paste0('pt', cova$Pseudotime)
cova$Pseudotime <- factor(x = cova$Pseudotime, levels =
                              c('pt1','pt2','pt3','pt4','pt5','pt6','pt7','pt0'))

# make sure everything matches
nrow(cova); ncol(log.norm.counts)
log.norm.counts <- log.norm.counts[ , colnames(log.norm.counts) %in% cova$RNAseqID ]
cova            <- cova[ cova$RNAseqID %in% colnames(log.norm.counts),  ]
nrow(cova) == ncol(log.norm.counts)
rownames( cova ) == colnames( log.norm.counts )

# CARNIVAL prior knowledge ------------------------------------------------

# We select interactions for human
AllInteractions = import_Omnipath_Interactions(select_organism = 9606)

# We transform to the format needed by CARNIVAL. We just keep signed and
# directed interactions
SignedDirectedInteractions <- filter(AllInteractions, is_directed == 1) %>%
    filter(is_stimulation == 1 | is_inhibition == 1)

InputCarnival <- bind_rows(
    (SignedDirectedInteractions %>%
         filter(is_stimulation == 1 & is_inhibition == 0) %>%
         transmute(source_genesymbol, interaction = 1, target_genesymbol)),
    (SignedDirectedInteractions %>%
         filter(is_stimulation == 0 & is_inhibition == 1) %>%
         transmute(source_genesymbol, interaction = -1, target_genesymbol))) %>%
    distinct()

# We have to be careful with the gene names with a "-". CPLEX gets crazy.
InputCarnival$source_genesymbol <-
    gsub("-","_",InputCarnival$source_genesymbol)
InputCarnival$target_genesymbol <-
    gsub("-","_",InputCarnival$target_genesymbol)

bad_int <- which(duplicated(paste(InputCarnival[,1],InputCarnival[,3])))

if ( length(bad_int) == 0 ) {
    InputCarnival <- InputCarnival
} else {
    InputCarnival = InputCarnival[-bad_int,]
}

# TF activities/Dorothea --------------------------------------------------

# Function to group Dorothea regulons.
# Input: A data frame containing Dorothea regulons, as stored in
# https://github.com/saezlab/ConservedFootprints/tree/master/data
# Output: Object of class regulon. See viper package.

df2regulon = function(df) {
    regulon = df %>%
        split(.$tf) %>%
        map(function(dat) {
            tf = dat %>% distinct(tf) %>% pull()
            targets = setNames(dat$mor, dat$target)
            likelihood = dat$likelihood
            list(tfmode =targets, likelihood = likelihood)
        })
    return(regulon)
}

# We read Dorothea Regulons for Human:
# https://github.com/saezlab/ConservedFootprints/tree/master/data/dorothea_benchmark/regulons
dorothea_regulon_human = read.csv("https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv",
                                  stringsAsFactors = F)

# We obtain the regulons based on interactions with confidence level A, B and C
regulon = dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C")) %>%
    df2regulon()
log.norm.counts_tfs = viper(eset = log.norm.counts, regulon = regulon, verbose = T, minsize = 15)
dim(log.norm.counts_tfs)

ctrl_ids <- which( cova$Pseudotime=="pt0" )

# Pseudotime 1 ------------------------------------------------------------

test_ids <- which( cova$Pseudotime=="pt1" )

# obtain t and p-values for pairwise comparisons
tmp <- viper::rowTtest( x = log.norm.counts_tfs[ , test_ids ], y = log.norm.counts_tfs[ , ctrl_ids ] )
doro <- data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2), p.value = signif(tmp$p.value, 3) )

# run CARNIVAL on t values from viper
target <- data.frame(t((doro$t)))
rownames(target) = "NES"
colnames(target) <- doro$Gene

CarnivalResults <- runCARNIVAL(netObj = InputCarnival,
                              measObj = target,
                              DOTfig = TRUE, dir_name = "CARNIVAL_Results_minsize_15",
                              threads = 10, #timelimit = 600,
                              solver = "cplex",
                              solverPath = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex")

write.table(CarnivalResults$weightedSIF,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_pt1_vs_pt0.csv", quote = F)
write.table(CarnivalResults$nodesAttributes,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_nodes_pt1_vs_pt0.csv", quote = F)
saveRDS(CarnivalResults, file = "CARNIVAL_Results_minsize_15/carnival_results_pt1_vs_pt0.rds")
saveRDS(doro, file = "CARNIVAL_Results_minsize_15/doro_pt1_vs_pt0.rds")

#CarnivalResults$sifAll[[1]]

# Pseudotime 2 ------------------------------------------------------------

test_ids <- which( cova$Pseudotime=="pt2" )

# obtain t and p-values for pairwise comparisons
tmp <- viper::rowTtest( x = log.norm.counts_tfs[ , test_ids ], y = log.norm.counts_tfs[ , ctrl_ids ] )
doro <- data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2), p.value = signif(tmp$p.value, 3) )

# run CARNIVAL on t values from viper
target <- data.frame(t((doro$t)))
rownames(target) = "NES"
colnames(target) <- doro$Gene

CarnivalResults <- runCARNIVAL(netObj = InputCarnival,
                               measObj = target,
                               DOTfig = TRUE, dir_name = "CARNIVAL_Results_minsize_15",
                               threads = 10, #timelimit = 600,
                               solver = "cplex",
                               solverPath = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex")

write.table(CarnivalResults$weightedSIF,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_pt2_vs_pt0.csv", quote = F)
write.table(CarnivalResults$nodesAttributes,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_nodes_pt2_vs_pt0.csv", quote = F)
saveRDS(CarnivalResults, file = "CARNIVAL_Results_minsize_15/carnival_results_pt2_vs_pt0.rds")
saveRDS(doro, file = "CARNIVAL_Results_minsize_15/doro_pt2_vs_pt0.rds")

# Pseudotime 3 ------------------------------------------------------------

test_ids <- which( cova$Pseudotime=="pt3" )

# obtain t and p-values for pairwise comparisons
tmp <- viper::rowTtest( x = log.norm.counts_tfs[ , test_ids ], y = log.norm.counts_tfs[ , ctrl_ids ] )
doro <- data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2), p.value = signif(tmp$p.value, 3) )

# run CARNIVAL on t values from viper
target <- data.frame(t((doro$t)))
rownames(target) = "NES"
colnames(target) <- doro$Gene

CarnivalResults <- runCARNIVAL(netObj = InputCarnival,
                               measObj = target,
                               DOTfig = TRUE, dir_name = "CARNIVAL_Results_minsize_15",
                               threads = 10, #timelimit = 600,
                               solver = "cplex",
                               solverPath = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex")

write.table(CarnivalResults$weightedSIF,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_pt3_vs_pt0.csv", quote = F)
write.table(CarnivalResults$nodesAttributes,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_nodes_pt3_vs_pt0.csv", quote = F)
saveRDS(CarnivalResults, file = "CARNIVAL_Results_minsize_15/CARNIVAL_Results_minsize_15_pt3_vs_pt0.rds")
saveRDS(doro, file = "CARNIVAL_Results_minsize_15/doro_pt3_vs_pt0.rds")

# Pseudotime 4 ------------------------------------------------------------

test_ids <- which( cova$Pseudotime=="pt4" )

# obtain t and p-values for pairwise comparisons
tmp <- viper::rowTtest( x = log.norm.counts_tfs[ , test_ids ], y = log.norm.counts_tfs[ , ctrl_ids ] )
doro <- data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2), p.value = signif(tmp$p.value, 3) )

# run CARNIVAL on t values from viper
target <- data.frame(t((doro$t)))
rownames(target) = "NES"
colnames(target) <- doro$Gene

CarnivalResults <- runCARNIVAL(netObj = InputCarnival,
                               measObj = target,
                               DOTfig = TRUE, dir_name = "CARNIVAL_Results_minsize_15",
                               threads = 10, #timelimit = 600,
                               solver = "cplex",
                               solverPath = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex")

write.table(CarnivalResults$weightedSIF,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_pt4_vs_pt0.csv", quote = F)
write.table(CarnivalResults$nodesAttributes,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_nodes_pt4_vs_pt0.csv", quote = F)
saveRDS(CarnivalResults, file = "CARNIVAL_Results_minsize_15/CARNIVAL_Results_minsize_15_pt4_vs_pt0.rds")
saveRDS(doro, file = "CARNIVAL_Results_minsize_15/doro_pt4_vs_pt0.rds")

# Pseudotime 5 ------------------------------------------------------------

test_ids <- which( cova$Pseudotime=="pt5" )

# obtain t and p-values for pairwise comparisons
tmp <- viper::rowTtest( x = log.norm.counts_tfs[ , test_ids ], y = log.norm.counts_tfs[ , ctrl_ids ] )
doro <- data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2), p.value = signif(tmp$p.value, 3) )

# run CARNIVAL on t values from viper
target <- data.frame(t((doro$t)))
rownames(target) = "NES"
colnames(target) <- doro$Gene

CarnivalResults <- runCARNIVAL(netObj = InputCarnival,
                               measObj = target,
                               DOTfig = TRUE, dir_name = "CARNIVAL_Results_minsize_15",
                               threads = 10, #timelimit = 600,
                               solver = "cplex",
                               solverPath = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex")

write.table(CarnivalResults$weightedSIF,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_pt5_vs_pt0.csv", quote = F)
write.table(CarnivalResults$nodesAttributes,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_nodes_pt5_vs_pt0.csv", quote = F)
saveRDS(CarnivalResults, file = "CARNIVAL_Results_minsize_15/CARNIVAL_Results_minsize_15_pt5_vs_pt0.rds")
saveRDS(doro, file = "CARNIVAL_Results_minsize_15/doro_pt5_vs_pt0.rds")

# Pseudotime 6 ------------------------------------------------------------

test_ids <- which( cova$Pseudotime=="pt6" )

# obtain t and p-values for pairwise comparisons
tmp <- viper::rowTtest( x = log.norm.counts_tfs[ , test_ids ], y = log.norm.counts_tfs[ , ctrl_ids ] )
doro <- data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2), p.value = signif(tmp$p.value, 3) )

# run CARNIVAL on t values from viper
target <- data.frame(t((doro$t)))
rownames(target) = "NES"
colnames(target) <- doro$Gene

CarnivalResults <- runCARNIVAL(netObj = InputCarnival,
                               measObj = target,
                               DOTfig = TRUE, dir_name = "CARNIVAL_Results_minsize_15",
                               threads = 10, #timelimit = 600,
                               solver = "cplex",
                               solverPath = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex")

write.table(CarnivalResults$weightedSIF,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_pt6_vs_pt0.csv", quote = F)
write.table(CarnivalResults$nodesAttributes,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_nodes_pt6_vs_pt0.csv", quote = F)
saveRDS(CarnivalResults, file = "CARNIVAL_Results_minsize_15/CARNIVAL_Results_minsize_15_pt6_vs_pt0.rds")
saveRDS(doro, file = "CARNIVAL_Results_minsize_15/doro_pt6_vs_pt0.rds")

# Pseudotime 7 ------------------------------------------------------------

test_ids <- which( cova$Pseudotime=="pt7" )

# obtain t and p-values for pairwise comparisons
tmp <- viper::rowTtest( x = log.norm.counts_tfs[ , test_ids ], y = log.norm.counts_tfs[ , ctrl_ids ] )
doro <- data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2), p.value = signif(tmp$p.value, 3) )

# run CARNIVAL on t values from viper
target <- data.frame(t((doro$t)))
rownames(target) = "NES"
colnames(target) <- doro$Gene

CarnivalResults <- runCARNIVAL(netObj = InputCarnival,
                               measObj = target,
                               DOTfig = TRUE, dir_name = "CARNIVAL_Results_minsize_15",
                               threads = 10, #timelimit = 600,
                               solver = "cplex",
                               solverPath = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex")

write.table(CarnivalResults$weightedSIF,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_pt7_vs_pt0.csv", quote = F)
write.table(CarnivalResults$nodesAttributes,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_Results_minsize_15/network_nodes_pt7_vs_pt0.csv", quote = F)
saveRDS(CarnivalResults, file = "CARNIVAL_Results_minsize_15/CARNIVAL_Results_minsize_15_pt7_vs_pt0.rds")
saveRDS(doro, file = "CARNIVAL_Results_minsize_15/doro_pt7_vs_pt0.rds")

