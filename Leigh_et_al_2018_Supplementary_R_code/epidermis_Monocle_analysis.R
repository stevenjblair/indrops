##Pseudotime of epidermis with Monocle

#homeostatic epidermis
library(monocle)
library(Seurat)

load("intact.Rdata")

#add IDs to metadata in preparation of importCDS
merged@meta.data$stashed.id <- merged@ident

pdf('intact_timecourse_tsne.pdf')
tsne.plot(merged, do.label = T)
dev.off()

WE_filtered <- SubsetData(merged, ident.remove = c(1,2,4,5,6,8,11,12))

pdf('intact_WE_timecourse_tsne.pdf')
tsne.plot(WE_filtered, do.label = T)
dev.off()

WE_filtered <- importCDS(WE_filtered)

WE_filtered <- detectGenes(WE_filtered, min_expr = .1)
expressed_genes <- row.names(subset(fData(WE_filtered),
                                    num_cells_expressed >= 3))

WE_filtered <- WE_filtered[expressed_genes,]

head(pData(WE_filtered))

WE_filtered <- estimateSizeFactors(WE_filtered)
WE_filtered <- detectGenes(WE_filtered, min_expr = .1)


head(fData(WE_filtered))
head(pData(WE_filtered))

pdf('intact_WE_variance.pdf')
plot_pc_variance_explained(WE_filtered, return_all = FALSE)
dev.off()

WE_filtered <- reduceDimension(WE_filtered,
                               max_components = 2,
                               norm_method = 'log',
                               num_dim = 8,
                               cores = 4,
                               residualModelFormulaStr = "~num_genes_expressed",
                               reduction_method = 'tSNE',
                               verbose = TRUE)

WE_filtered <- clusterCells(WE_filtered, verbose = TRUE)

WE_filtered$tree.ident <- as.character(WE_filtered$tree.ident)

pdf('intact_WE_clusters.pdf')
plot_cell_clusters(WE_filtered)
dev.off()

pdf('intact_WE_numGenes.pdf')
plot_cell_clusters(WE_filtered, color_by = "num_genes_expressed")
dev.off()

pdf('intact_WE_byseurat.pdf')
plot_cell_clusters(WE_filtered, color_by = "stashed.id")
dev.off()


pdf('intact_WE_rho_delta.pdf')
plot_rho_delta(WE_filtered, rho_threshold = 2, delta_threshold = 10)
dev.off()

clustering_DEG_genes <- differentialGeneTest(WE_filtered,
                                             fullModelFormulaStr = '~stashed.id',
                                             cores = 4)

my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:200]
WE_filtered <- setOrderingFilter(WE_filtered, ordering_genes = my_ordering_genes)
WE_filtered <- reduceDimension(WE_filtered,
                               method = 'DDRTree',
                               fullModelFormulaStr = '~stashed.id', verbose = F, scaling = T, maxIter = 100, norm_method = 'log', max_components = 6, param.gamma = 20)

WE_filtered <- orderCells(WE_filtered)

pdf('intact_WE_pseudotime.pdf')
plot_cell_trajectory(WE_filtered, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_bySeuratCluster.pdf')
plot_cell_trajectory(WE_filtered, color_by = "stashed.id", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_byCluster.pdf')
plot_cell_trajectory(WE_filtered, color_by = "Cluster", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_bySeuratCluster_facet.pdf')
plot_cell_trajectory(WE_filtered, color_by = "stashed.id") +
  facet_wrap(~stashed.id, nrow = 1)
dev.off()

save.image("082018_intact_WE_postDEG.Rdata")

###

load("082018_intact_WE_postDEG.Rdata")

# create vector of no's
my_vector <- rep('no', nrow(pData(WE_filtered)))

# change status to yes if the cell was in cluster 10
my_vector[pData(WE_filtered)$stashed.id == 10] <- rep('yes', sum(pData(WE_filtered)$stashed.id == 10))

# add vector to phenoData
pData(WE_filtered)$test <- my_vector
head(pData(WE_filtered))

clustering_DEG_genes <- differentialGeneTest(WE_filtered,
                                             fullModelFormulaStr = '~test',
                                             cores = 4)

my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:150]
WE_filtered <- setOrderingFilter(WE_filtered, ordering_genes = my_ordering_genes)
WE_filtered <- reduceDimension(WE_filtered,
                               method = 'DDRTree',
                               fullModelFormulaStr = '~stashed.id', verbose = F, scaling = T, maxIter = 100, norm_method = 'log', max_components = 15, param.gamma = 20)

WE_filtered <- orderCells(WE_filtered, reverse = TRUE)

pdf('intact_WE_pseudotime.pdf')
plot_cell_trajectory(WE_filtered, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_bySeuratCluster.pdf')
plot_cell_trajectory(WE_filtered, color_by = "stashed.id", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_byCluster.pdf')
plot_cell_trajectory(WE_filtered, color_by = "Cluster", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_bySeuratCluster_facet.pdf')
plot_cell_trajectory(WE_filtered, color_by = "stashed.id", show_branch_points = FALSE) +
  facet_wrap(~stashed.id, nrow = 1)
dev.off()

###

save.image("intact_WE_postDEG.Rdata")

rm(list=setdiff(ls(), "WE_filtered")) 

save.image("intact_WE_slim.Rdata")

png('intact_WE_facet.png', width = (7/4)*3, height = 2, units = 'in', res = 300)
plot_cell_trajectory(WE_filtered, color_by = "stashed.id", show_branch_points = FALSE) +
  facet_wrap(~stashed.id, nrow = 1)+
  theme(legend.position="none")
dev.off()

png('intact_WE_pseudotime.png', width = 2.1, height = 1.77, units = 'in', res = 300)
plot_cell_trajectory(WE_filtered, color_by = "Pseudotime", show_branch_points = FALSE)+
  theme(legend.position="none")
dev.off()

###Regenerating epidermis


library(monocle)
library(Seurat)

#load medium-bud blastema data (medium-bud and 23dpa are equivalent)
load("medium_bud_blastema.RData")

pdf('timecourse_tsne.pdf')
tsne.plot(combined, do.label = T)
dev.off()

combined <- importCDS(combined)

pData(combined)$"day1" <- pData(combined)$protocol == "d1" 

combined <- detectGenes(combined, min_expr = .1)
expressed_genes <- row.names(subset(fData(combined),
                                    num_cells_expressed >= 3))

#only want expressed genes and the day 1 subset of cells
combined_day1 <- combined[expressed_genes, pData(combined)$"day1"]

#only want WE populations
pData(combined_day1)$"WE" <- pData(combined_day1)$tree.ident == "1" |
  pData(combined_day1)$tree.ident == "2" |
  pData(combined_day1)$tree.ident == "3" |
  pData(combined_day1)$tree.ident == "4"

combined_day1_WE <- combined_day1[expressed_genes, pData(combined_day1)$"WE"]

head(pData(combined_day1_WE))

combined_day1_WE <- estimateSizeFactors(combined_day1_WE)
combined_day1_WE <- detectGenes(combined_day1_WE, min_expr = .1)


head(fData(combined_day1_WE))
head(pData(combined_day1_WE))

pdf('23dpa_1_WE_variance.pdf')
plot_pc_variance_explained(combined_day1_WE, return_all = FALSE)
dev.off()

combined_day1_WE <- reduceDimension(combined_day1_WE,
                                    max_components = 2,
                                    norm_method = 'log',
                                    num_dim = 8,
                                    cores = 4,
                                    residualModelFormulaStr = "~num_genes_expressed",
                                    reduction_method = 'tSNE',
                                    verbose = TRUE)

combined_day1_WE <- clusterCells(combined_day1_WE, verbose = TRUE)

#plot clusters
pdf('23dpa_day1_WE_clusters.pdf')
plot_cell_clusters(combined_day1_WE)
dev.off()

pdf('23dpa_day1_WE_numGenes.pdf')
plot_cell_clusters(combined_day1_WE, color_by = "num_genes_expressed")
dev.off()

pdf('23dpa_day1_WE_byseurat.pdf')
plot_cell_clusters(combined_day1_WE, color_by = "tree.ident")
dev.off()


pdf('WE_rho_delta.pdf')
plot_rho_delta(combined_day1_WE, rho_threshold = 2, delta_threshold = 10)
dev.off()

clustering_DEG_genes <- differentialGeneTest(combined_day1_WE,
                                             fullModelFormulaStr = '~tree.ident',
                                             cores = 4)

my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:150]
combined_day1_WE <- setOrderingFilter(combined_day1_WE, ordering_genes = my_ordering_genes)
combined_day1_WE <- reduceDimension(combined_day1_WE,
                                    method = 'DDRTree',
                                    fullModelFormulaStr = '~tree.ident', verbose = F, scaling = T, norm_method = 'log', max_components =8, param.gamma = 20)

combined_day1_WE$tree.ident <- as.character(combined_day1_WE$tree.ident)

combined_day1_WE <- orderCells(combined_day1_WE)

#plot pseudotime graphs
pdf('23dpa_day1_WE_pseudotime.pdf')
plot_cell_trajectory(combined_day1_WE, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

pdf('23dpa_day1_WE_bySeuratCluster.pdf')
plot_cell_trajectory(combined_day1_WE, color_by = "tree.ident", show_branch_points = FALSE)
dev.off()

pdf('23dpa_day1_WE_byState.pdf')
plot_cell_trajectory(combined_day1_WE, color_by = "Cluster", show_branch_points = FALSE)
dev.off()

pdf('23dpa_day1_WE_bySeuratCluster_facet.pdf')
plot_cell_trajectory(combined_day1_WE, color_by = "tree.ident", show_branch_points = FALSE) +
  facet_wrap(~tree.ident, nrow = 1)
dev.off()

save.image("23dpa_day1_WE_postDEG.Rdata")

rm(list=setdiff(ls(), "combined_day1_WE")) 

save.image("23dpa_day1_WE_slim.Rdata")

png('23dpa_day1_WE_facet.png', width = (7/5)*3, height = 2, units = 'in', res = 300)
plot_cell_trajectory(combined_day1_WE, color_by = "tree.ident", show_branch_points = FALSE) +
  facet_wrap(~tree.ident, nrow = 1)+
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(-6,-3,0,2))
dev.off()

png('23dpa_day1_WE_pseudotime.png', width = (7/5)*1.15, height = 1.70, units = 'in', res = 300)
plot_cell_trajectory(combined_day1_WE, color_by = "Pseudotime", show_branch_points = FALSE) +
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(-6,-3,0,2))
dev.off()

png('23dpa_day1_WE_pseudotime_K17.png', width = 7, height = 2, units = 'in', res = 300)
genes_to_plot <- c('c1083200_g2_i3^sp|A1L595|K1C17_BOVIN^Filament')
plot_genes_in_pseudotime(combined_day1_WE[genes_to_plot,], color_by = "tree.ident")
dev.off()



