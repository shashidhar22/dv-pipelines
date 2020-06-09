library(monocle3)
library(Seurat)
library(BiocSingular)
library(scran)
library(scater)
library(batchelor)
library(patchwork)
library(tidyverse)
library(DropletUtils)
set.seed(12357)
sampler <- function(sce) {
  sce_path <- "/Volumes/home/fuser/nextflow/MMRFINE/Counts"
  print(paste(sce_path, sce, "outs/filtered_feature_bc_matrix/", sep = "/"))
  sce <- read10xCounts(paste(sce_path, sce, "outs/filtered_feature_bc_matrix/", sep="/"))
  sce <- sce[,sample(ncol(sce), 500, replace = TRUE)]
  return(sce)
}
#Load samples and subsample the, combine them and cluster them without correction
samples <- c( "201687_6B_0", "202823_6P_0", "211146_18P_0", "256339_22P_1", "293138_60P_0", "325592_60P_1")
sce_list <- samples %>% map(sampler) %>% map(logNormCounts)
gene_var <- sce_list %>% map(modelGeneVar)
gene_comb <- do.call(combineVar, gene_var)
gene_hvgs <- gene_comb$bio > 0
sce_comb <- do.call(cbind, sce_list)
sce_comb <- runPCA(sce_comb, subset_row=gene_hvgs, BSPARAM=BiocSingular::RandomParam())
sce_comb <- runUMAP(sce_comb, dimred="PCA")
sce_uncorrected <- plotUMAP(sce_comb, colour_by=as.factor("Batch")) + ggtitle("Scater uncorrected")

#Run multibatch correction
sce_comb <- do.call(multiBatchNorm, sce_list)
sce_comb <- do.call(cbind, sce_list)
sce_comb <- runPCA(sce_comb, subset_row=gene_hvgs, BSPARAM=BiocSingular::RandomParam())
sce_comb <- runUMAP(sce_comb, dimred="PCA")
sce_multibatch <- plotUMAP(sce_comb, colour_by=as.factor("Batch")) + ggtitle("Scater multibatch corrected")

#Run FastMNN
sce_comb <- do.call(cbind, sce_list)
sce_comb$batch <- 1
sce_comb <- fastMNN(sce_comb, batch=sce_comb$Batch, subset.row = gene_hvgs, BSPARAM=BiocSingular::RandomParam(deferred = TRUE))
dim(reducedDim(sce_comb, "corrected"))
sce_comb <- runUMAP(sce_comb, dimred="corrected")
sce_fastMNN <- plotUMAP(sce_comb, colour_by="batch") + ggtitle("fastMNN")

(sce_fastMNN | sce_multibatch | sce_uncorrected)

cell_types <- c("S1_Bcells" = 1, "S1_CytotoxicTCells" = 2, "S1_HelperTCells" = 3, "S1_MemoryTCells" = 4,
                "S1_Monocytes" = 5, "S1_NaiveCytotoxicTCells" = 6, "S1_NaiveTCells" = 7, "S1_NKCells" = 8,
                "S1_Progenitor" = 9, "S1_RegulatoryTCells" = 10)

# tenX aggregate
sce_aggr <- read10xCounts("Counts/Aggregate_mapped_normalized/outs/filtered_feature_bc_matrix/")
sample_sheet <- c("1" = "S1_Bcells", "2" = "S1_CytotoxicTCells", "3" = "S1_HelperTCells", "4" = "S1_MemoryTCells", "5" = "S1_Monocytes", "6" = "S1_NaiveCytotoxicTCells",
                  "7" = "S1_NaiveTCells", "8" = "S1_NKCells", "9" = "S1_Progenitor", "10" = "S1_RegulatoryTCells")
colData(sce_aggr) <- colData(sce_aggr) %>% as_tibble() %>% mutate(Sample=str_extract(Barcode,pattern = "\\d+$")) %>% mutate(Sample= recode(Sample, !!!sample_sheet)) %>% DataFrame()
sce_aggr <- logNormCounts(sce_aggr)
gene_var <- modelGeneVar(sce_aggr)
variable_genes <- getTopHVGs(gene_var, n=1000)
sce_aggr <- runPCA(sce_aggr, subset_row=variable_genes)
pcs <- reducedDim(sce_aggr)
choices <- getClusteredPCs(sce_aggr)
reducedDim(sce_aggr, "PCA_chosen") <- reducedDim(sce_aggr, "PCA")[,1:metadata(choices)$chosen]
sce_aggr <- runUMAP(sce_aggr, dimred="PCA")
sce_aggrplot <- plotUMAP(sce_aggr, colour_by="Sample") + ggtitle("10X Aggregated")
#((sce_fastMNN + sce_multibatch +sce_moncle) / (sce_uncorrected + p2  + sce_aggrplot))

#Seurat integrate
seu_sampler <- function(sce) {
  sce <- readRDS(sce)
  sample <- sce$Sample[1]
  #sce <- sce[,sample(ncol(sce), 500, replace = TRUE)]
  sce <- SubsetData(sce, max.cells.per.ident = 1500, random.seed = 12357)
  sce <- NormalizeData(sce)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 1000)
  sce$Barcode <- sce$Barcode %>% str_replace("1", paste('1', as.character(cell_types[sce$Sample[1]]), sep = '_')) %>% as_vector()
  return(sce)
}
sce_list <- list.files("Preprocess/CDS", pattern = "*_filtered_seurat_cds.rds", full.names = TRUE) %>% map(seu_sampler) 
immune.anchors <- FindIntegrationAnchors(object.list = sce_list, dims = 1:50)
integrated <- IntegrateData(anchorset = immune.anchors, dims = 1:50)
DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:50)
sce_seuratStd <- DimPlot(integrated, reduction = "umap", group.by = "Sample", label = FALSE, repel = TRUE) + ggtitle("Suerat standard protocol")

#Run sctranform
for (i in 1:length(sce_list)) {
  sce_list[[i]] <- SCTransform(sce_list[[i]], verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = sce_list, nfeatures = 1000)
sce_list <- PrepSCTIntegration(object.list = sce_list, anchor.features = features, 
                                    verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = sce_list, normalization.method = "SCT", 
                                           anchor.features = features, verbose = FALSE)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:50)
sce_seuratSCT <- DimPlot(integrated, reduction = "umap", group.by = "Sample", label = TRUE, repel = TRUE) + ggtitle("Seurat SCTransform")

#Monocle combine
mon_sampler <- function(sce) {
  sce <- readRDS(sce)
  sample <- sce$Sample[1]
  sce <- sce[,sample(ncol(sce), 1500, replace = TRUE)]
  return(sce)
}
sce_list <- list.files("Preprocess/CDS", pattern = "*_filtered_monocle3_cds.rds", full.names = TRUE) %>% map(readRDS)
sce_list <- combine_cds(sce_list)
sce_list <- preprocess_cds(sce_list, method = "PCA", norm_method = "log")
sce_list <- align_cds(cds = sce_list, num_dim = 50, alignment_group = "Sample", preprocess_method = "PCA")
sce_list <- reduce_dimension(sce_list, reduction_method = "UMAP", preprocess_method = "PCA")
sce_moncle <- plot_cells(sce_list, color_cells_by="Sample", label_cell_groups=FALSE) + ggtitle("Monocle integration")
((sce_fastMNN + sce_multibatch +sce_moncle) / (sce_uncorrected + sce_seuratSCT + sce_seuratStd ))

