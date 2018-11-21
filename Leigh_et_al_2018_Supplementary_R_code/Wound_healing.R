#Wound healing samples
#This section analyses three samples collected during wound healing (3 days post amputation) 

#First, load the required packages. 
library(Seurat)
library(dplyr)

#load in data
data = read.table('wound_healing_and_medium_bud.repGene', header=T, row.names=1, sep='\t')
#pull out samples N4, N5, and N6 which are the wound healing stage limb samples
N4_N5_N6 = data[,grep('^N[456]', colnames(data))]

#Remove data matrix with extra samples
rm(data)
#Create Seurat object and make sparse
seurat_3dpa = CreateSeuratObject(N4_N5_N6, project = '3dpa', min.cells = 8, min.genes = 200)
seurat_3dpa = MakeSparse(seurat_3dpa)

#While we did some filtering above, we need to perform further quality control to ensure that the cells we are working with aren't apoptotic #or have a dearth of genes. First, we need to identify the mitochondrial genes present in this matrix. The axolotl mitochondrial genome can #be found here: https://www.ncbi.nlm.nih.gov/nuccore/AJ584639. Remember that the genes are written as protein names when greping for #mitochondrial genes. 

#find mitochonrial genes in matrix. The protein name should be used and changed for each gene within the mitochondrial genome.
grep(pattern = "*CYB_*", x = rownames(x = seurat_3dpa.intact@data), value = TRUE)
#list of all mitochondrial genes in this wound healing matrix
mito.genes.3dpa <- c("c786641_g1_i1^sp|Q9B205|CYB_CAICR", "c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1043008_g2_i1^sp|Q9B205|CYB_CAICR", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c1057599_g1_i2^sp|P00018|CYC_DRONO^Cytochrome_CBB3", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c1053715_g3_i1^sp|O03539|COX1_NOTPE", "c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c289614_g1_i1^sp|P05503|COX1_RAT", "c959712_g1_i1^sp|P00416|COX3_MOUSE", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1025234_g1_i1^sp|O63796|NU2M_ANACA", "c1068681_g4_i4^sp|Q9ZZM3|NU5M_SALSA")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.3dpa <- Matrix::colSums(seurat_3dpa@raw.data[mito.genes.3dpa, ])/Matrix::colSums(seurat_3dpa@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_3dpa <- AddMetaData(object = seurat_3dpa, metadata = percent.mito.3dpa, col.name = "percent.mito")

#Now perform quality control on matrix by filtering out cells with high percent mitochondrial RNA and low and high number of genes. We filter #out cells that by visualize inspection appear to have relatively high mitochondrial RNA content or high or low number of genes. These #numbers can be modified to be more or less inclusive. 

#visualize number of genes, unique molecular identifiers (UMI), and percent mitochondrial RNA
VlnPlot(object = seurat_3dpa, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#filter out cells
seurat_3dpa <- FilterCells(object = seurat_3dpa, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(5000, 0.15))

#normalize data
seurat_3dpa <- NormalizeData(object = seurat_3dpa, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable genes
seurat_3dpa <- FindVariableGenes(object = seurat_3dpa, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#scale data and regress out nUMI and percent.mito
seurat_3dpa <- ScaleData(object = seurat_3dpa, vars.to.regress = c("nUMI", "percent.mito"))


#Next, we perform linear dimensional reduction and visualize the results in a few different ways. 
seurat_3dpa <- RunPCA(object = seurat_3dpa, pc.genes = seurat_3dpa@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#visualize results
VizPCA(object = seurat_3dpa, pcs.use = 1:2)
PCAPlot(object = seurat_3dpa, dim.1 = 1, dim.2 = 2)
seurat_3dpa <- ProjectPCA(object = seurat_3dpa, do.print = FALSE)
PCHeatmap(object = seurat_3dpa, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = seurat_3dpa, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = seurat_3dpa)

#plot standard deviations to chose PCs to use in downstream analysis, here we chose 19
seurat_3dpa <- FindClusters(object = seurat_3dpa, reduction.type = "pca", dims.use = 1:19, resolution = 1.5, print.output = 0, save.SNN = TRUE)

#Now we can identify cell populations within the homeostatic limb, vizualize the resulting populations using tSNE, and subsequently find markers that define these different populations 

#find clusters using first 18 PCs
seurat_inDrops3.intact <- FindClusters(object = seurat_inDrops3.intact, reduction.type = "pca", dims.use = 1:18, resolution = 1.5, print.output = 0, save.SNN = TRUE)

#run non-linear dimensional reduction
seurat_3dpa <- RunTSNE(object = seurat_3dpa, dims.use = 1:19, do.fast = TRUE)

# Build a phylogenetic tree to see how cells are related while simultaneously renaming and reordering cluster names according to their #position on the tree. This will be important to determine when deciding whether similar populations should be merged. 
seurat_3dpa <- BuildClusterTree(seurat_3dpa, do.reorder=TRUE, reorder.numeric=TRUE)

#visualize tSNE 
set.seed(5)
TSNEPlot(object = seurat_3dpa)
#visulize tSNE based on sample to determine how similar the two samples are to one another
TSNEPlot(object = seurat_3dpa, group.by = "orig.ident")

#assess nodes
node.scores <- AssessNodes(seurat_3dpa)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores
node.scores


#merge first 5 nodes
#select nodes to merge
nodes.merge=node.scores[1:5,]
nodes.to.merge <- sort(x = nodes.merge$node)

#create a new Seurat object in which we will merge our selected nodes
merged <- seurat_3dpa
#merge nodes
for (n in nodes.to.merge) {merged <- MergeNode(object = merged, node.use = n)}


#re-visualize the tSNE after we have merged the non-distinct nodes
set.seed(5)
TSNEPlot(merged, do.label = TRUE)


#determine differentially expressed genes for each population

#find markers for each population
all.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#write DE results to table for inspection
write.table(all.markers, 'wound.healing.markers.txt', sep = '\t')

#save Rdata
save.image('wound_healing.Rdata')

#end session
q()