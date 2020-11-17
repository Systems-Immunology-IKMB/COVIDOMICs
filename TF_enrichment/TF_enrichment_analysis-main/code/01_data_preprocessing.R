#
# @author:	Axel Künstner
# @email:	axel.kuenstner@uni-luebeck.de
# @date:	2020-07-08
#

library(tidyverse)
library(skimr)
# Get data ----------------------------------------------------------------

cova <- read.delim("../data/Sequenced_sample_metadata_020720.csv", sep = ",")

norm.counts <- read.delim("../data/COVID_control_normalized_counts_070720.txt")

# Reformat covariates -----------------------------------------------------

head(cova)

# remove samples without RNAseqID
cova <- cova %>%
#    tibble %>%
    dplyr::filter( RNAseqID != "" ) %>%
    dplyr::mutate( Patient = gsub(pattern = "#", replacement = "P", Patient) ) %>%
    dplyr::mutate( PatientID = factor(PatientID) ) %>%
    dplyr::mutate( Remission = factor(Remission) ) %>%
    dplyr::mutate( Gender = factor(Gender) ) %>%
    dplyr::mutate( Patient = factor(Patient) ) %>%
    dplyr::mutate( Disease.traj = factor(Disease.traj) ) %>%
    dplyr::mutate( Pseudotime = factor(Pseudotime) )

rownames(cova) <- cova$RNAseqID

cova <- cova[ order( rownames( cova ) ), ]

# check data
skimr::skim_without_charts(cova)

table( cova$Disease.traj )
table( cova$Disease.traj, cova$Pseudotime )

saveRDS(object = cova, file = "../data/covariates.RDS")

# Reformat Expression data ------------------------------------------------

# add gene symbols to expression data
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
t2g  <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
head(t2g)

norm.counts <- merge( norm.counts, t2g, by.x = "row.names", by.y = "ensembl_gene_id" )
norm.counts$Row.names <- NULL

# sort by expression AND
# remove duplicated genes, keep copy with highest expression
idx <- order( rowSums( norm.counts[ 1: (ncol(norm.counts)-1)],  ), decreasing = TRUE)
norm.counts <- norm.counts[idx, ]
idx <- which( !duplicated(norm.counts$external_gene_name) )
norm.counts <- norm.counts[ idx, ]
# just a check
norm.counts[ norm.counts$external_gene_name == "FAM95B1", ]

rownames(norm.counts) <- norm.counts$external_gene_name
norm.counts$external_gene_name <- NULL

norm.counts <- norm.counts[ , order( colnames(norm.counts) ) ]

saveRDS(object = norm.counts, file = "../data/RNAseq_cnts_norm.RDS")

