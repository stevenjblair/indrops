#Pseudotime analysis of blastema-resident cells

#load required packages
library(URD)
library(Seurat)
library(dplyr)
#load medium-bud data
load('medium_bud_blastema.RData')

#pull out all non-immune blastema cells
dpa23_FDB <- SubsetData(combined.merged, ident.use=c('5','6','10', '13', '14','15', '16','17'), subset.name = 'day', accept.value = 'd1')
#remove medium-bud Seurat object
rm(combined.merged)

#add time point metadata to Seurat object
dpa23_FDB@meta.data$time <- "medium-bud"

#load in early-bud data
load('early_bud_blastema.Rdata')

#pull out all non-immune blastema cells
dpa14_blastema <- SubsetData(merged, ident.use= c('3', '4', '5','6','11'))
#add time point metadata to Seurat object
dpa14_blastema@meta.data$time <- "early-bud"

#remove excess data
rm(seurat_14dpa, percent.mito.14dpa, mito.genes.14dpa, merged)

#load in wound healing data
load('wound_healing.Rdata')
#pull out all non-immune blastema cells
dpa3_blastema <- SubsetData(merged, ident.use = c('2','21','24','25'))

#add time point metadata to Seurat object
dpa3_blastema@meta.data$time <- 'wound_healing'

#remove excess data
rm(seurat_3dpa, mito.genes.3dpa, percent.mito.3dpa, node.scores, merged, all.markers, nodes.merge, nodes.to.merge, n)

#combine above Seurat object so we can make a raw data file and metadata file with only the cells we want to put into URD
dpa14.23 <- MergeSeurat(dpa14_blastema, dpa23_FDB)
combined <- MergeSeurat(dpa14.23,dpa3_blastema)

#pull out raw data
raw.data <- as.matrix(combined@raw.data)

#pull out metadata
meta.data <- as.matrix(combined@meta.data)

#remove excess Seurat objects
rm(dpa14.23, dpa3_blastema, dpa14_blastema, dpa23_FDB)

#create URD object with all of the cells we selected above
blastema <- createURD(count.data = raw.data, meta = meta.data, min.cells=3, min.counts=3)

#now remove raw data and Seurat objects
rm(raw.data, meta.data, combined)

#add time point info to URD stage slot
blastema@group.ids$stage <- as.character(blastema@meta[rownames(blastema@group.ids),"time"])

#find variable genes for each time point
stages <- sort(unique(blastema@group.ids$stage))
var.3dpa <- findVariableGenes(blastema, cells.fit=cellsInCluster(blastema, "stage", 'wound_healing'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.14dpa <- findVariableGenes(blastema, cells.fit=cellsInCluster(blastema, "stage", 'early-bud'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.23dpa <- findVariableGenes(blastema, cells.fit=cellsInCluster(blastema, "stage", 'medium-bud'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)

#combine into one 
var.genes <- sort(unique(unlist(var.3dpa, var.14dpa, var.23dpa)))

#add to URD object
blastema@var.genes <- var.genes

#calculate PCA
blastema <- calcPCA(blastema, mp.factor = 2)

pcSDPlot(blastema)

# Calculate tSNE
set.seed(19)
blastema <- calcTsne(object = blastema)

#visualize tSNE by time point
plotDim(blastema, "time", plot.title = "tSNE: DPA")

#create an URD object with only medium-bud cells
blastema.23dpa <- urdSubset(blastema, cells.keep=cellsInCluster(blastema, "stage", "medium-bud"))

#use variable genes that were calculated only for medium-bud
blastema.23dpa@var.genes <- var.23dpa

#calculate PCA and tSNE
blastema.23dpa <- calcPCA(blastema.23dpa, mp.factor = 1.5)
pcSDPlot(blastema.23dpa)
set.seed(20)
blastema.23dpa <- calcTsne(blastema.23dpa)

#perform graph based clustering
blastema.23dpa <- graphClustering(blastema.23dpa, num.nn = c(30, 40, 50, 70), method = 'Louvain', do.jaccard=T)
blastema.23dpa <- graphClustering(blastema.23dpa, num.nn = c(50, 60, 70, 80, 100), method = 'Infomap', do.jaccard=T)
clusterings <- c(paste0('Louvain-',c(30, 40, 50, 70)), paste0('Infomap-', c(50, 60, 70, 80, 100)))
for (c in clusterings) {plot(plotDim(blastema.23dpa, c, legend=F))}

#we chose to go forward with Infomap-100
clusters <- unique(blastema.23dpa@group.ids$'Infomap-100')

#find markers for these populations to get an idea of what they are
pr.markers <- lapply(clusters, function(c) markersAUCPR(blastema.23dpa, clust.1 = c, clustering = 'Infomap-100', genes.use= blastema.23dpa@var.genes))

#you can export each table like below, or if R studio look at top markers in viewer
#note that cluster number and number in pr.markers[[]] are not necessarily equal!
#cluster number can be found in headings of the marker lists
write.table(pr.markers[[1]], 'clus1.markers.txt', sep = '\t')

#looking at cluster genes to see if any cells we don't want to include may have made it through prior clustering
head(pr.markers[[3]], 20)

#we will not use clusters 7 (Myeloid cells), 10 (Wound epidermis), 13 (Erythrocytes) so these aren't included in trajectory analysis
#get cell names for the cells we will use in URD
dpa23_good_cells <- cellsInCluster(blastema.23dpa, clustering = 'Infomap-100', cluster = c('1','2','3','4','5','6','8','9','11','12'))

#let's clean up the other two time points since it's clear that WE, erythrocytes, etc. may be contaminating the clusters
#create URD object with just early-bud sample
blastema.14dpa <- urdSubset(blastema, cells.keep=cellsInCluster(blastema, "stage", "early-bud"))

#use variable genes found for early-bud
blastema.14dpa@var.genes <- var.14dpa

#calculate PCA and tSNE
blastema.14dpa <- calcPCA(blastema.14dpa, mp.factor = 1.5)
pcSDPlot(blastema.14dpa)
set.seed(20)
blastema.14dpa <- calcTsne(blastema.14dpa)

#perform graph-based clustering
blastema.14dpa <- graphClustering(blastema.14dpa, num.nn = c(30, 40, 50, 70), method = 'Louvain', do.jaccard=T)
blastema.14dpa <- graphClustering(blastema.14dpa, num.nn = c(50, 60, 70, 80, 100), method = 'Infomap', do.jaccard=T)
clusterings <- c(paste0('Louvain-',c(30, 40, 50, 70)), paste0('Infomap-', c(50, 60, 70, 80, 100)))
for (c in clusterings) {plot(plotDim(blastema.14dpa, c, legend=F))}

#iwe chose to go forward with Infomap-50 based clustering
clusters <- unique(blastema.14dpa@group.ids$'Infomap-50')

#calculate marker genes
pr.markers_14dpa <- lapply(clusters, function(c) markersAUCPR(blastema.14dpa, clust.1 = c, clustering = 'Infomap-50', genes.use= blastema.14dpa@var.genes))

#inspect genes that define subsets, remembering that cluster number is found in column names
#either or export to table or inspect in R/Rstudio
write.table(pr.markers_14dpa[[1]], '081218.clus1.txt', sep = '\t')
#looking at cluster genes
head(pr.markers_14dpa[[10]])

#need to remove clusters 1 (WE), 7 (WE), 11(Erythrocytes), 13(Immune) which would be 217 cells that likely wont participate in this
#get names of cells we will use in URD
dpa14_good_cells <- cellsInCluster(blastema.14dpa, clustering = 'Infomap-50', cluster = c('2','3','4','5','6','8','9','10','12','14','15','16')) 

#finally, let's clean the wound healing cells
#create a subsetted object of cells from wound healing
blastema.3dpa <- urdSubset(blastema, cells.keep=cellsInCluster(blastema, "stage", "wound_healing"))

#use wound healing variable genes
blastema.3dpa@var.genes <- var.14dpa

#calculate PCA and tSNE
blastema.3dpa <- calcPCA(blastema.3dpa, mp.factor = 1.5)

pcSDPlot(blastema.3dpa)
set.seed(20)
blastema.3dpa <- calcTsne(blastema.3dpa)

#graph-based clustering
blastema.3dpa <- graphClustering(blastema.3dpa, num.nn = c(30, 40, 50, 70), method = 'Louvain', do.jaccard=T)
blastema.3dpa <- graphClustering(blastema.3dpa, num.nn = c(50, 60, 70, 80, 100), method = 'Infomap', do.jaccard=T)
clusterings <- c(paste0('Louvain-',c(30, 40, 50, 70)), paste0('Infomap-', c(50, 60, 70, 80, 100)))
for (c in clusterings) {plot(plotDim(blastema.3dpa, c, legend=F))}


#we chose to go forward with Infomap-60
clusters <- unique(blastema.3dpa@group.ids$'Infomap-60')

#find marker genes
pr.markers_3dpa <- lapply(clusters, function(c) markersAUCPR(blastema.3dpa, clust.1 = c, clustering = 'Infomap-60', genes.use= blastema.3dpa@var.genes))

#inspect marker genes to determine what to remove
head(pr.markers_3dpa[[8]])

#we will remove 3 (WE) and 12 (erythrocytes)
dpa3_good_cells <- cellsInCluster(blastema.3dpa, clustering = 'Infomap-60', cluster = c('1','2','4','5','6','7','8','9','10','11','13'))

#collect all cells to use in URD
paste(dpa23_good_cells, dpa14_good_cells) -> cells_to_use
paste(cells_to_use, dpa3_good_cells) -> cells_for_urd


#make an URD object with the cells we want to use
cleaned <- urdSubset(blastema, cells.keep= c(dpa3_good_cells, dpa23_good_cells, dpa14_good_cells))

#re-calculate variable genes for these
stages <- sort(unique(cleaned@group.ids$stage))
var.3dpa <- findVariableGenes(cleaned, cells.fit=cellsInCluster(cleaned, "stage", 'wound_healing'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.14dpa <- findVariableGenes(cleaned, cells.fit=cellsInCluster(cleaned, "stage", 'early-bud'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.23dpa <- findVariableGenes(cleaned, cells.fit=cellsInCluster(cleaned, "stage", 'medium-bud'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.genes <- sort(unique(unlist(var.3dpa, var.14dpa, var.23dpa)))
cleaned@var.genes <- var.genes

#calculate PCA
cleaned <- calcPCA(cleaned, mp.factor = 2)
pcSDPlot(cleaned)

#calculate tSNE, see Figure 7a
set.seed(19)
cleaned <- calcTsne(object = cleaned)
plotDim(cleaned, "time", plot.title = "tSNE by time point", legend = F, point.size = 2)

#calculate diffision map, allowing destinty to determine sigma (this value was determined to be 19.538, which we used)
cleaned <- calcDM(cleaned, knn = 54, sigma = 19.538)

#visualize dim arrays
plotDimArray(cleaned, reduction.use = "dm", dims.to.plot = 1:16, outer.title = "Diffusion Map (Sigma 19.538, 54 NNs): DPA", label="stage", plot.title="", legend=F)

#tsne with transitions
plotDim(cleaned, "time", transitions.plot = 10000, plot.title="DPA (with transitions)")

#use cells from wound healing as root
root.cells <- cellsInCluster(cleaned, "stage", "wound_healing")

#run 'flood' simulations
cleaned.floods <- floodPseudotime(cleaned, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

#process simulations into a pseudotime
cleaned <- floodPseudotimeProcess(cleaned, cleaned.floods, floods.name="pseudotime")

#check for adequate number of simulations (should reach asymptote)
pseudotimePlotStabilityOverall(cleaned)

#visualize tSNE with pseudotime overlaid
plotDim(cleaned, "pseudotime")

#plot pseudtime at each time point, Figure 7b
plotDists(cleaned, "pseudotime", "time", plot.title="Pseudotime by time point", legend = F)

#create URD object with just medium-bud cells
cleaned.23dpa <- urdSubset(cleaned, cells.keep=cellsInCluster(cleaned, "stage", "medium-bud"))

#use medium-bud variable genes
cleaned.23dpa@var.genes <- var.23dpa

# Calculate PCA and tSNE
cleaned.23dpa <- calcPCA(cleaned.23dpa, mp.factor = 1.5)
pcSDPlot(cleaned.23dpa)
set.seed(20)
cleaned.23dpa <- calcTsne(cleaned.23dpa)

#perform graph-based clustering
cleaned.23dpa <- graphClustering(cleaned.23dpa, num.nn = c(30, 40, 50, 70), method = 'Louvain', do.jaccard=T)
cleaned.23dpa <- graphClustering(cleaned.23dpa, num.nn = c(50, 60, 70, 80, 100), method = 'Infomap', do.jaccard=T)
clusterings <- c(paste0('Louvain-',c(30, 40, 50, 70)), paste0('Infomap-', c(50, 60, 70, 80, 100)))
for (c in clusterings) {plot(plotDim(cleaned.23dpa, c, legend=T))}

#tSNE plot Figure 7c
plotDim(cleaned.23dpa, 'Infomap-100', plot.title = "Medium-bud blastema populations", legend = F, point.size = 2)

#we chose to go forward with Infomap-100
clusters <- unique(cleaned.23dpa@group.ids$'Infomap-100')

#find marker genes, see Supplementary Data files for marker lists
pr.markers_23dpa_cleaned <- lapply(clusters, function(c) markersAUCPR(cleaned.23dpa, clust.1 = c, clustering = 'Infomap-100', genes.use= cleaned.23dpa@var.genes))
#inspect markers in R
head(pr.markers_23dpa_cleaned[[3]], 20)

#distal blastema population has some WE markers, let's remove WE cells
plotDot(cleaned.23dpa, genes = c('c1069858_g1_i1^sp|Q90X25|HXA13_CHICK^HoxA13_N','c1070920_g1_i4^sp|Q2VL56|PAX9_SAGOE^PAX^Tm2','c1020768_g1_i2^sp|Q9H2S6|TNMD_HUMAN^BRICHOS^Tm1','c1083312_g1_i2^sp|P70390|SHOX2_MOUSE','c1081900_g1_i4^sp|P25815|S100P_HUMAN','c1091168_g1_i2^sp|Q66S13|NATT4_THANI^DUF946^sigP'), clustering = 'Infomap-100')
distal.score <- apply(cleaned.23dpa@logupx.data[c('c1069858_g1_i1^sp|Q90X25|HXA13_CHICK^HoxA13_N','c1070920_g1_i4^sp|Q2VL56|PAX9_SAGOE^PAX^Tm2','c1020768_g1_i2^sp|Q9H2S6|TNMD_HUMAN^BRICHOS^Tm1','c1083312_g1_i2^sp|P70390|SHOX2_MOUSE'), cellsInCluster(cleaned.23dpa, 'Infomap-100', '1')], 2, sum.of.logs)
new.distal <- names(which(distal.score > 0))
remove.distal <- names(which(distal.score <= 0))

#add names for all pops
i100.n <- length(unique(cleaned.23dpa@group.ids$'Infomap-100'))
i100.cluster.assignments <- data.frame(clusters = 1:i100.n, name = rep(NA, i100.n), tip = rep(NA, i100.n), row.names = 1:i100.n)
i100.cluster.assignments['1', 'name'] <- "Distal Blastema"
i100.cluster.assignments['2', 'name'] <- "FAPs"
i100.cluster.assignments['3', 'name'] <- "Synovial Fibroblasts"
i100.cluster.assignments['4', 'name'] <- "Cartilage"
i100.cluster.assignments['5', 'name'] <- "Osteoblast-like"
i100.cluster.assignments['6', 'name'] <- "Joint"
i100.cluster.assignments['7', 'name'] <- "Schwann" 
i100.cluster.assignments['8', 'name'] <- "Endothelial"
i100.cluster.assignments['9', 'name'] <- "Myogenic"
i100.cluster.assignments['10', 'name'] <- "Pericytes" 


cluster.assignments <- i100.cluster.assignments

cluster.assignments <- cluster.assignments[!is.na(cluster.assignments$name), ]
cluster.assignments$cluster.new <- 1:nrow(cluster.assignments)


cleaned@group.ids$'23dpa-Infomap-100' <- NA
cleaned@group.ids[rownames(cleaned.23dpa@group.ids), '23dpa-Infomap-100'] <- cleaned.23dpa@group.ids$'Infomap-100'
cleaned@group.ids$'23dpa-Cluster' <- NA
cleaned.23dpa@group.ids$clusters.23dpa.name <- NA
cleaned.23dpa@group.ids$clusters.23dpa.num <- NA

#so now for cluster 1 I'll just give it the cleaned new.distal 
cleaned.23dpa@group.ids[new.distal, "clusters.23dpa.name"] <- cluster.assignments[1, "name"] 
cleaned.23dpa@group.ids[new.distal, "clusters.23dpa.num"] <- as.character(cluster.assignments[1, "cluster.new"])

cells_2 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '2') 
cleaned.23dpa@group.ids[cells_2, "clusters.23dpa.name"] <- cluster.assignments[2, "name"] 
cleaned.23dpa@group.ids[cells_2, "clusters.23dpa.num"] <- as.character(cluster.assignments[2, "cluster.new"])

cells_3 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '3') 
cleaned.23dpa@group.ids[cells_3, "clusters.23dpa.name"] <- cluster.assignments[3, "name"] 
cleaned.23dpa@group.ids[cells_3, "clusters.23dpa.num"] <- as.character(cluster.assignments[3, "cluster.new"])

cells_4 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '4') 
cleaned.23dpa@group.ids[cells_4, "clusters.23dpa.name"] <- cluster.assignments[4, "name"] 
cleaned.23dpa@group.ids[cells_4, "clusters.23dpa.num"] <- as.character(cluster.assignments[4, "cluster.new"])

cells_5 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '5') 
cleaned.23dpa@group.ids[cells_5, "clusters.23dpa.name"] <- cluster.assignments[5, "name"] 
cleaned.23dpa@group.ids[cells_5, "clusters.23dpa.num"] <- as.character(cluster.assignments[5, "cluster.new"])

cells_6 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '6') 
cleaned.23dpa@group.ids[cells_6, "clusters.23dpa.name"] <- cluster.assignments[6, "name"] 
cleaned.23dpa@group.ids[cells_6, "clusters.23dpa.num"] <- as.character(cluster.assignments[6, "cluster.new"])

cells_7 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '7') 
cleaned.23dpa@group.ids[cells_7, "clusters.23dpa.name"] <- cluster.assignments[7, "name"] 
cleaned.23dpa@group.ids[cells_7, "clusters.23dpa.num"] <- as.character(cluster.assignments[7, "cluster.new"])

cells_8 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '8') 
cleaned.23dpa@group.ids[cells_8, "clusters.23dpa.name"] <- cluster.assignments[8, "name"] 
cleaned.23dpa@group.ids[cells_8, "clusters.23dpa.num"] <- as.character(cluster.assignments[8, "cluster.new"])

cells_9 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '9') 
cleaned.23dpa@group.ids[cells_9, "clusters.23dpa.name"] <- cluster.assignments[9, "name"] 
cleaned.23dpa@group.ids[cells_9, "clusters.23dpa.num"] <- as.character(cluster.assignments[9, "cluster.new"])

cells_10 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '10') 
cleaned.23dpa@group.ids[cells_10, "clusters.23dpa.name"] <- cluster.assignments[10, "name"] 
cleaned.23dpa@group.ids[cells_10, "clusters.23dpa.num"] <- as.character(cluster.assignments[10, "cluster.new"])


cleaned@group.ids$'23dpa-Infomap-100' <- NA
cleaned@group.ids[rownames(cleaned.23dpa@group.ids), '23dpa-Infomap-100'] <- cleaned.23dpa@group.ids$'Infomap-100'
cleaned@group.ids[rownames(cleaned.23dpa@group.ids), '23dpa-Cluster'] <- cleaned.23dpa@group.ids$clusters.23dpa.name
cleaned@group.ids$'23dpa-Cluster-Num'<- NA
cleaned@group.ids[rownames(cleaned.23dpa@group.ids), '23dpa-Cluster-Num'] <- cleaned.23dpa@group.ids$clusters.23dpa.num


cleaned@group.ids[rownames(cleaned.23dpa@group.ids), "tip.clusters"] <- cleaned.23dpa@group.ids$clusters.23dpa.num

#determine potential for terminal populations
potential <- clusterTipPotential(cleaned, 'pseudotime', 'tip.clusters', name.store = 'tip.potential')
potential 

#tips = 3, 4, 5, 6, 7, 8, 9, 10, so everything but "distal blastema" and "FAPs"

only.tips <- cellsInCluster(cleaned, clustering= '23dpa-Cluster-Num', cluster = c('3','4','5','6','7','8','9','10'))
tips <- urdSubset(cleaned, cells.keep= only.tips)
cleaned@group.ids[rownames(tips@group.ids), "real.tip.clusters"] <- tips@group.ids$'23dpa-Cluster-Num'

#determine parameters of the logistic used to bias the transition probabilities
cleaned.ptlogistic <- pseudotimeDetermineLogistic(cleaned, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)

#bias the transition matrix acording to pseudotime
cleaned.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(cleaned, "pseudotime", logistic.params=cleaned.ptlogistic))

#simulate the biased random walks from each tip
cleaned.walks <- simulateRandomWalksFromTips(cleaned, tip.group.id="real.tip.clusters", root.cells=root.cells, transition.matrix = cleaned.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)

#process the biased random walks into visitation frequencies
cleaned <- processRandomWalksFromTips(cleaned, cleaned.walks, verbose = F)

#color only tip clusters on tSNE
plotDim(cleaned, "real.tip.clusters", plot.title="Cells in each tip")

#load tip clusters into tree
cleaned.tree <- loadTipCells(cleaned, "real.tip.clusters")

#build tree
cleaned.tree <- buildTree(cleaned.tree, pseudotime = "pseudotime", tips.use=c('3','4','5','6','7','8','9','10'), divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)

#rename clusters
cleaned.tree <- nameSegments(cleaned.tree, segments=c('3','4','5','6','7','8','9','10'), segment.names = c("Synovial Fibros", "Cartilage", "Osteoblast-like", "Joint","Schwann", "Endothelial", "Myogenic", "Pericyte"), short.names = c("Synovial Fibros", "Cartilage",  "Osteoblast-like", "Joint","Schwann", "Endothelial", "Myogenic", "Pericyte"))

#plot tree with time point info overlaid
plotTree(cleaned.tree, "stage", title="DPA")

#plot tree with medium-bud blastema colors overlaid
plotTree(cleaned.tree, 'tip.clusters', title = '23dpa_clusters', cell.alpha = 0.5, cell.size = 2, tree.alpha = 0.5, tree.size = .25)

joint.markers <- aucprTestAlongTree(cleaned.tree, pseudotime="pseudotime", tips='Joint', log.effect.size=0.4, auc.factor = 1.25, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use=genes.use, only.return.global=F, must.beat.sibs=0.6, report.debug=T)

synovial.markers <- aucprTestAlongTree(cleaned.tree, pseudotime="pseudotime", tips='Synovial Fibros', log.effect.size=0.4, auc.factor = 1.25, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use=genes.use, only.return.global=F, must.beat.sibs=0.6, report.debug=T)

#visuzalize population specific markers
#osteoblast-like

p1<- plotTree(cleaned.tree, 'c862122_g2_i1^sp|Q6DJ00|OSTCN_XENTR^sp|P40147|OSTCN_XENLA^Gla^sigP', title = 'OSTCN')
p2 <- plotTree(cleaned.tree, 'c1034953_g1_i2^sp|A1YQ92|ODAM_MACMU', title = 'ODAM')

#example of graph in Supplementary Fig7
plot_grid(p1,p2)

#can continue to overlay markers like shown above

#save
save.image('blastema.URD.RData')

#quit
q()