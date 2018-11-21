# Homeostatic limb
#This section analyses the homeostatic limb, which includes two samples. 

#First, load the required packages. 
library(Seurat)
library(dplyr)


#Now, load the cell by gene matrix that includes the two homeostatic limb samples. This cell by gene matrix has three other samples, that we #will not use right now. So we need to select out only samples 1 and 2. We will use intact to as shorthand for homeostatic limbs. 

#load in data
inDrops3.data = read.table('intact_and_contralateral.repGene', header = T, row.names = 1, sep = '\t')
#pull out samples 1 and 2, which are the intact limb samples
inDrops3.intact = inDrops3.data[,grep('^S[12]_', colnames(inDrops3.data))]

#Remove data matrix with extra samples
rm(inDrops3.data) 
#Create Seurat object and make sparse
seurat_inDrops3.intact = CreateSeuratObject(inDrops3.intact, project = 'inDrops3.intact', min.cells = 8, min.genes = 200)
seurat_inDrops3.intact = MakeSparse(seurat_inDrops3.intact)

#While we did some filtering above, we need to perform further quality control to ensure that the cells we are working with aren't apoptotic #or have a dearth of genes. First, we need to identify the mitochondrial genes present in this matrix. The axolotl mitochondrial genome can #be found here: https://www.ncbi.nlm.nih.gov/nuccore/AJ584639. Remember that the genes are written as protein names when greping for #mitochondrial genes. 

#find mitochonrial genes in matrix. The protein name should be used and changed for each gene within the mitochondrial genome.
grep(pattern = "*CYB_*", x = rownames(x = seurat_inDrops3.intact@data), value = TRUE)
#list of all mitochondrial genes in this intact matrix
mito.genes.intact <- c("c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1025234_g1_i1^sp|O63796|NU2M_ANACA", "c1068681_g4_i1^sp|Q9ZZM3|NU5M_SALSA^sp|P82013|VDAC2_MELGA^Porin_3", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.intact <- Matrix::colSums(seurat_inDrops3.intact@raw.data[mito.genes.intact, ])/Matrix::colSums(seurat_inDrops3.intact@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_inDrops3.intact <- AddMetaData(object = seurat_inDrops3.intact, metadata = percent.mito.intact, col.name = "percent.mito")

#Now perform quality control on matrix by filtering out cells with high percent mitochondrial RNA and low and high number of genes. We filter #out cells that by visualize inspection appear to have relatively high mitochondrial RNA content or high or low number of genes. These #numbers can be modified to be more or less inclusive. 

#visualize number of genes, unique molecular identifiers (UMI), and percent mitochondrial RNA
VlnPlot(object = seurat_inDrops3.intact, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#filter out cells
seurat_inDrops3.intact <- FilterCells(object = seurat_inDrops3.intact, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(4000, 0.125))

#normalize data
seurat_inDrops3.intact <- NormalizeData(seurat_inDrops3.intact, normalization.method= "LogNormalize", scale.factor= 10000)

#find variable genes
seurat_inDrops3.intact <- FindVariableGenes(object = seurat_inDrops3.intact, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)

#scale data and regress out nUMI and percent.mito
seurat_inDrops3.intact <- ScaleData(seurat_inDrops3.intact, vars.to.regress = c('nUMI', 'percent.mito'))


#Next, we perform linear dimensional reduction and visualize the results in a few different ways. 

seurat_inDrops3.intact <- RunPCA(object = seurat_inDrops3.intact, pc.genes = seurat_inDrops3.intact@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#visualize results
PrintPCA(object = seurat_inDrops3.intact, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = seurat_inDrops3.intact, pcs.use = 1:2)
PCAPlot(object = seurat_inDrops3.intact, dim.1 = 1, dim.2 = 2)
seurat_inDrops3.intact <- ProjectPCA(object = seurat_inDrops3.intact, do.print = FALSE)
PCHeatmap(object = seurat_inDrops3.intact, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = seurat_inDrops3.intact, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = seurat_inDrops3.intact, pc.use = 13:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

#plot standard deviations to chose PCs to use in downstream analysis, here we chose 18
PCElbowPlot(object = seurat_inDrops3.intact)


#Now we can identify cell populations within the homeostatic limb, vizualize the resulting populations using tSNE, and subsequently find markers that define these different populations 

#find clusters using first 18 PCs
seurat_inDrops3.intact <- FindClusters(object = seurat_inDrops3.intact, reduction.type = "pca", dims.use = 1:18, resolution = 1.5, print.output = 0, save.SNN = TRUE)

#run non-linear dimensional reduction
seurat_inDrops3.intact <- RunTSNE(object = seurat_inDrops3.intact, dims.use = 1:18, do.fast = TRUE)

# Build a phylogenetic tree to see how cells are related while simultaneously renaming and reordering cluster names according to their #position on the tree. This will be important to determine when deciding whether similar populations should be merged. 
seurat_inDrops3.intact <- BuildClusterTree(seurat_inDrops3.intact, do.reorder=TRUE, reorder.numeric=TRUE)

#visualize tSNE 
set.seed(5)
TSNEPlot(object = seurat_inDrops3.intact, do.label = T)

#visulize tSNE based on sample to determine how similar the two samples are to one another
TSNEPlot(object = seurat_inDrops3.intact, group.by = 'orig.ident')

#assess nodes
node.scores <- AssessNodes(seurat_inDrops3.intact)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores
node.scores


#merge first 7 nodes
#select nodes to merge
nodes.merge <- node.scores[1:7, ]
nodes.to.merge <- sort(x = nodes.merge$node)

#create a new Seurat object in which we will merge our selected nodes
merged <- seurat_inDrops3.intact
#merge nodes
for (n in nodes.to.merge) {merged <- MergeNode(object = merged, node.use = n)}


#re-visualize the tSNE after we have merged the non-distinct nodes
set.seed(5)
TSNEPlot(merged, do.label = TRUE)


#determine differentially expressed genes for each population

#find markers for each population
all.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#write DE results to table for inspection
write.table(all.markers, 'intact.only.markers.txt', sep = '\t')

#save Rdata
save.image('intact.Rdata')

#end session
q()
