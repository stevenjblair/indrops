#Medium-bud stage blastema samples
#This section analyses six samples collected during medium bud stage (23 days post amputation) 

#First, load the required packages. 
library(Seurat)
library(dplyr)

#load in first three medium-bud samples
data.d1 = read.table("early_and_medium_bud.repGene", header=T, row.names=1, sep='\t')
#sample S1, S2, and S4 are medium-bud samples collected on day 1
S1_S2_S4 = data.d1[,grep("^S[124]_",colnames(data.d1))]

#create Seurat object and make sparse
seurat_S1_S2_S4 = CreateSeuratObject(raw.data = S1_S2_S4, project = "23dpa_d1", min.cells = 8, min.genes = 200)
seurat_S1_S2_S4 <- MakeSparse(seurat_S1_S2_S4) 

#list mito genes medium-bud day 1 matrix
mito.genes.i4.d1 <- c("c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1043008_g2_i1^sp|Q9B205|CYB_CAICR", "c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c786641_g1_i1^sp|Q9B205|CYB_CAICR", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1057599_g1_i1^sp|P00018|CYC_DRONO^Cytochrome_CBB3", "c1127119_g1_i1^sp|P81280|CYC_ALLMI", "c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c959712_g1_i1^sp|P00416|COX3_MOUSE", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1060846_g2_i1^sp|P03921|NU5M_MOUSE", "c1068681_g4_i1^sp|Q9ZZM3|NU5M_SALSA^sp|P82013|VDAC2_MELGA^Porin_3")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.i4.d1 <- Matrix::colSums(seurat_S1_S2_S4@raw.data[mito.genes.i4.d1, ])/Matrix::colSums(seurat_S1_S2_S4@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_S1_S2_S4 <- AddMetaData(object = seurat_S1_S2_S4, metadata = percent.mito.i4.d1, col.name = "percent.mito") 

#remove raw data
rm(data.d1, S1_S2_S4)

#load in second set of three medium-bud samples collected on day 2
data.d2 = read.table('wound_healing_and_medium_bud.repGene', header=T, row.names=1, sep='\t')
#samples N1, N2, and N3 are medium-bud stage blastema samples collected on day 2
N1_N2_N3 = data.d2[,grep("^N[123]",colnames(data.d2))]

#create seurat object and make sparse
seurat_N1_N2_N3 = CreateSeuratObject(raw.data = N1_N2_N3, project = "23dpa_d2", min.cells = 8, min.genes = 200)
seurat_N1_N2_N3 <- MakeSparse(seurat_N1_N2_N3)

#list mito genes from day 2 matrix
mito.genes.i4.d2 <- c("c786641_g1_i1^sp|Q9B205|CYB_CAICR", "c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1043008_g2_i1^sp|Q9B205|CYB_CAICR", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c289614_g1_i1^sp|P05503|COX1_RAT", "c959712_g1_i1^sp|P00416|COX3_MOUSE", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1060846_g2_i1^sp|P03921|NU5M_MOUSE", "c1068681_g4_i4^sp|Q9ZZM3|NU5M_SALSA")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.i4.d2 <- Matrix::colSums(seurat_N1_N2_N3@raw.data[mito.genes.i4.d2, ])/Matrix::colSums(seurat_N1_N2_N3@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_N1_N2_N3 <- AddMetaData(object = seurat_N1_N2_N3, metadata = percent.mito.i4.d2, col.name = "percent.mito")

#save each sample info in a metadata sample slot
seurat_N1_N2_N3 <- StashIdent(seurat_N1_N2_N3, save.name = 'sample')
seurat_S1_S2_S4 <- StashIdent(seurat_S1_S2_S4, save.name = 'sample')

#custom filters
seurat_S1_S2_S4 <- FilterCells(object = seurat_S1_S2_S4, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(6000, 0.1))
seurat_N1_N2_N3 <- FilterCells(object = seurat_N1_N2_N3, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(6000, 0.1))

#normalize data
seurat_N1_N2_N3 <- NormalizeData(object = seurat_N1_N2_N3, normalization.method = "LogNormalize", scale.factor= 10000)
seurat_S1_S2_S4 <- NormalizeData(object = seurat_S1_S2_S4, normalization.method = "LogNormalize", scale.factor= 10000)

#find variable genes
seurat_N1_N2_N3 <- FindVariableGenes(object = seurat_N1_N2_N3, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
seurat_S1_S2_S4 <- FindVariableGenes(object = seurat_S1_S2_S4, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#scale data and regress out nUMI and percent.mito
seurat_N1_N2_N3 <- ScaleData(object = seurat_N1_N2_N3, vars.to.regress = c('nUMI', 'percent.mito'))
seurat_S1_S2_S4 <- ScaleData(object = seurat_S1_S2_S4, vars.to.regress = c('nUMI', 'percent.mito'))

#add day collected information into a metadata column
seurat_S1_S2_S4@meta.data$day <- "d1"
seurat_N1_N2_N3@meta.data$day <- "d2"

#find highly variable genes
hvg.d1 <- rownames(x = head(x = seurat_S1_S2_S4@hvg.info, n = 2000))
hvg.d2 <- rownames(x = head(x = seurat_N1_N2_N3@hvg.info, n = 2000))
hvg.union <- union(x = hvg.d1, y = hvg.d2)

#run canonical correlation analysis and save resulting matrix as a Seurat object called combined
combined <- RunCCA(object = seurat_S1_S2_S4, object2= seurat_N1_N2_N3, genes.use = hvg.union, num.cc = 30)

#determine cca dimensions to use in downstream analysis
DimHeatmap(object = combined, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
DimHeatmap(object = combined, reduction.type = "cca", cells.use = 500, dim.use = 10:19, do.balanced = TRUE)
DimHeatmap(object = combined, reduction.type = "cca", cells.use = 500, dim.use = 20:29, do.balanced = TRUE)

#we chose first 25 dimensions
combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "day", dims.align = 1:25)

#perform integrated analysis on all cells
combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:25, do.fast = T)
combined <- FindClusters(combined, reduction.type = "cca.aligned", resolution = 1, dims.use = 1:25)

# Build a phylogenetic tree to see how cells are related while simultaneously renaming and reordering cluster names according to their #position on the tree. This will be important to determine when deciding whether similar populations should be merged. 
combined <- BuildClusterTree(combined, do.reorder=TRUE, reorder.numeric=TRUE)

#assess nodes
node.scores <- AssessNodes(combined)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores
node.scores

#merge first 2 nodes
nodes.merge <- node.scores[1:2, ] 
nodes.to.merge <- sort(x = nodes.merge$node)
combined.merged <- combined 
for (n in nodes.to.merge) {combined.merged <- MergeNode(object = combined.merged, node.use = n) }

#visualize tSNE by sample, day collected, and assigned cluster 
TSNEPlot(combined.merged, do.return = T, pt.size = 0.5, group.by = "sample")
TSNEPlot(combined.merged, do.return = T, pt.size = 0.5, group.by = "day")
TSNEPlot(combined.merged, do.label = T, do.return = T, pt.size = 0.5)

#find DE genes for each population
all.markers.merged <- FindAllMarkers(object = combined.merged, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, logfc.threshold = 0.35, max.cells.per.ident = 2000)

#write DE results to table for inspection
write.table(all.markers.merged, 'all.markers.medium.bud.txt', sep = '\t') 

save.image('medium_bud_blastema.RData')

q()