#
# @author:	Axel Künstner / Hauke Busch
# @email:	axel.kuenstner@uni-luebeck.de
# @date:	2020-08-13
#
# Dorothea/CARNIVAL workflow was done according to
# https://github.com/saezlab/NetworkModeling_course_2020
#

# following https://github.com/saezlab/NetworkModeling_course_2020

library(tidyverse)
library(ggplot2)
library(dnet)
library(igraph)
library(EnhancedVolcano)
library(pheatmap)
library(CINNA)

# Get TF activity and p-value data ----------------------------------------------------------------

data_directory <- "CARNIVAL_Results_minsize_15"

timepoints     <- paste("pt",1:7,sep="")
dorothea_files <- paste(data_directory,"/","doro_",timepoints,"_vs_pt0.rds",sep="")


dorothea_results<- list()

for(i in 1:length(timepoints)){
	dorothea_results[[timepoints[i]]] <- readRDS(dorothea_files[i])
}

# Prepare data matrices for Heatmap and Volcano plot -----------------------

dorothea_t <- cbind(dorothea_results$pt1$t, # t statistics
				  	dorothea_results$pt2$t,
				  	dorothea_results$pt3$t,
				  	dorothea_results$pt4$t,
				  	dorothea_results$pt5$t,
				  	dorothea_results$pt6$t,
				  	dorothea_results$pt7$t
				  	)

rownames(dorothea_t) <- rownames(dorothea_results$pt1) # Should all be the same. but maybe include a check on this
colnames(dorothea_t) <- timepoints

dorothea_p <- cbind(dorothea_results$pt1$p, # p values
				  	dorothea_results$pt2$p,
				  	dorothea_results$pt3$p,
				  	dorothea_results$pt4$p,
				  	dorothea_results$pt5$p,
				  	dorothea_results$pt6$p,
				  	dorothea_results$pt7$p
				  	)

rownames(dorothea_p) <- rownames(dorothea_results$pt1) # Should all be the same. but maybe include a check on this
colnames(dorothea_p) <- timepoints

# Filter data and plot as heatmap

# P-value cutoff - either plot * or nothing
max_pvalue     <- 0.05
dorothea_p <- ifelse(dorothea_p < max_pvalue,"*"," ")

# With percent of the most varying TFs to display
iqr_threshold <- 0.50
dorothea_iqr  <- apply(dorothea_t,1,IQR)
t_threshold   <- quantile(myiqr,iqr_threshold)

# Maximal t-statitic to plot
max_tstatistic <- 8
dorothea_t[dorothea_t >  dthr] <-  dthr
dorothea_t[dorothea_t < -dthr] <- -dthr

pheatmap(dorothea_t[names(dorothea_iqr)[dorothea_iqr > t_threshold],],
         cellwidth=10,
         cellheight=5,
         border_color=NA,
         cluster_cols=F,
         fontsize=5,
         display_number=dorothea_p[names(dorothea_iqr)[dorothea_iqr > t_threshold],]
         )



# Get the output from the Carinval simulation ---------------------------------------------------------

carnival_files <- paste(data_directory,"/","carnical_results_",timepoints,"_vs_pt0.rds",sep="")

# A list to hold the network information at the different time points
carnival_networks <- list()

# Loop of the carnival output
k <- 1
for(ff in myfiles){

	tmplist <- list()

	# Get data at Pseudotime
	CarnivalResults <- readRDS(file =ff)
	net   			<- data.frame(CarnivalResults$weightedSI,stringsAsFactors = F)
	nodes 			<- data.frame(CarnivalResults$nodesAttributes, stringsAsFactors = F) %>%
    							  mutate(act_sign = sign(as.numeric(AvgAct))
    							  )
    rownames(nodes) <- nodes$Node

	# Create network with weight cutoff 50
	net.adj 		<- net %>%
    					dplyr::mutate( Weight = as.numeric(Weight) ) %>%
    					dplyr::filter( Weight >= 50 ) %>%
    					dplyr::filter( Node1 != "Perturbation" ) %>%
    					dplyr::mutate( edge.color = "black" )
    net.adj$edge.color[ net.adj$Sign == -1 ] <- "red"

    # Encode the edges
    net.adj$eid 		<-	apply(net.adj,1,function(x){paste(x[1],x[3],sep="|")})
    rownames(net.adj) 	<- net.adj$eid
    net.adj 			<- select(net.adj, -c(eid))

    # Encode Node Attributes and store as igraph network
	cgr 				<- net.adj %>%
    						dplyr::filter( Weight >= 50 ) %>%
    						dplyr::select( c("Node1", "Node2","Weight") ) %>%
    						igraph::graph_from_data_frame(directed = TRUE, vertices = NULL)
    V(cgr)$color 		<- ifelse(nodes[V(cgr)$name,"act_sign"] >0,"#E6E6FA","#FFE4E1")

    # Store objects temporarily
	tmplist[["net"]] 		 <- net
	tmplist[["nodes"]]		 <- nodes
	tmplist[["net.adj"]]	 <- net.adj
	tmplist[["cgr"]] 		 <- cgr

	# Store in carnival_networks
	carnival_networks[[timepoints[k]]] <- tmplist
	k <- k+1
}

# Calculate and plot the degree centrality   ---------------------------------------

my_centrality <- "Degree Centrality"
# Not all nodes are present in all timepoints, so first get the union of all node names
all_nodes <- NULL
for( i in timepoints){
	all_nodes <- c(carnival_networks(mytime[[i]]$net.adj),all_nodes)
}
all_nodes <- sort(unique(all_nodes))

# Create a data structure holding the degree centrality
# Initialized with 0
centrality 				<- as.data.frame(mat.or.vec(length(all_nodes),length(timepoints))
names(centrality) 		<- timepoints
rownames(centrality) 	<- all_nodes

# Loop over all timepoints
for(pts in timepoints){
	cgr  <- carnival_networks[[pts]]$cgr
	# Extract the largest connected component
	cgr_c <- dNetInduce(g=cgr,
						nodes_query=V(cgr)$name,
						knn=0,
						remove.loops=F,
						largest.comp=T
						)
	# Which centrality measures are proper?
	pr_cent	<-	proper_centralities(cgr_c)
	idx 	<- which(pr_cent == "Degree Centrality")

	deg_centrality <- unlist(calculate_centralities(cgr_c,
								include = pr_cent[idx]
								)
							)
	# Make the result a named vector
	names(deg_centrality) <- V(cgr_c)$name

	# And store the results
	centrality[names(deg_centrality),pts] <- deg_centrality
}


# Plot the centrality as a heatmap

# With min. IQR to display
iqr_threshold <- 1

# Exclude nodes that have low overall centrality
sum_cutoff <- 10

# Display cutoff in the heatmap
display_cutoff <- 12

centality_iqr 	<- apply(centrality,IQR)
centrality_sum 	<- apply(centrality,1,sum)

# Select the nodes
idx <- which(centrality_sum > sum_cutoff &
			centality_iqr > iqr_threshold)

# Display cutoff
centrality[centrality > display_cutoff] <- display_cutoff

# Plot the heatmap
pheatmap(centrality[idx,],
		cluster_col=F,
		cellwidth=10,
		cellheight=5,
		fontsize=5,
		border_color=NA
		)



# Plot nearest neighbor network around MAPK1 and MAPK3 for pseudotime 3 -------------------------

# Define timepoint
pts <- timepoint[3]

# Get network object
cgr  	<- 	 carnival_networks[[pts]]$cgr
net.adj <-   carnival_networks[[pts]]$net.adj

# Color the edges
for(i in 1:nrow(net.adj)){
	 	 E(cgr)$color[i] <- net.adj[i,"edge.color"]
	 }


# Construct a network around MAPK3 and MAPK1 with 2 nearest neighbors
# Retain only the largest component
cgr_c <- dNetInduce(g=cgr,
                    nodes_query=c("MAPK3","MAPK1"),
                    knn=2,
                    remove.loops=F,
                    largest.comp=T)

# Scale node size according to degree
vsize <- 0.7*degree(cgr_c)+4

# Output filename
fname <- paste("SubNetwork",pts,".pdf",sep="_")

# Plot Network
pdf(file = fname, height = 16, width = 12 )
visNet(cgr_c,
		vertex.shape="circle",
		newpage=F,
        vertex.label.cex=1,
        vertex.color = V(cgr_c)$color,
        vertex.size=vsize,
        edge.color=E(cgr_c)$color,
        glayout =  layout.fruchterman.reingold(cgr_c,niter = 5000),
        edge.arrow.size=1,
        edge.curved=0.1,
        vertex.frame.color="grey50"
        )
dev.off()


# Volcano plot of TFs in patients with or without remission -------------------------------------------
# Here we use the regulon data defined in 02_CARNIVAL_preprocessing

data_directory <- "../data/"

# Load the metadata from all three cohorts
covariate_file <- paste(data_directory,"/All_cohorts_covariates.rds",sep="")
covariates_all_cohorts <- readRDS(covariate_file)

# Load the normalized count data for all three cohorts
count_file <- paste(data_directory,"/All_cohorts_normalized_counts.rds",sep="")
cnts.norm <- readRDS(covariate_file)

# Infer the TF activity using viper
norm.counts_tfs = viper(eset = cnts.norm, regulon = regulon, verbose = T, minsize = 15)


# Compare Remission versus no remission in pseudotime 2 ---------------------------------------

# Test IDs
test_ids <- which( covariates_all_cohorts$Pseudotime=="pt2" &
                   covariates_all_cohorts$Remission=="NR" &
                   covariates_all_cohorts$Site=="Nederland")
# Control IDs
ctrl_ids <- which( 	covariates_all_cohorts$Pseudotime=="pt2" &
 					covariates_all_cohorts$Remission=="R" &
 					covariates_all_cohorts$Site=="Nederland")

# obtain t and p-values for pairwise comparisons
tmp <- viper::rowTtest( x = norm.counts_tfs[ , test_ids ], y = norm.counts_tfs[ , ctrl_ids ] )
dorothea_nr_vs_r <- data.frame(Gene = rownames(tmp$p.value),
								t = round(tmp$statistic, 2),
								p.value = signif(tmp$p.value, 2) )

# Volcano plot
EnhancedVolcano(dorothea_nr_vs_r,
    	lab = rownames(dorothea_nr_vs_r),
    	x = 't',
    	y = 'p.value',
    	xlim = c(-4, 4),
        ylim = c(0, 4),
        legendPosition = 'right',
        pCutoff = 0.1,
        FCcutoff = 1,
        pointSize = 1.0,
        drawConnectors = TRUE,
        widthConnectors = 0.2,
        colConnectors = 'grey30',
        title = 'No remission vs remission',
         subtitle = ' ',
         xlab = "t-statistics"
         )

