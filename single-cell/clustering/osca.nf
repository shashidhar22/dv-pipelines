#!/usr/bin/env nextflow
// Create input channels
cds_path = Channel.fromPath(params.input.cds_path + "/*_sce.rds")
adownsample = Channel.from(params.osca.abinitio.nsample)
anfeatures = Channel.from(params.osca.abinitio.nfeatures)
acenters = Channel.from(params.osca.abinitio.centers)
rmethod = Channel.from(params.osca.abinitio.reduction_method)
cmethod = Channel.from(params.osca.abinitio.clustering_method)
study_id = Channel.from(params.input.study_name)
// Get log counts
process OSCA_LOGCOUNT{
  echo false
  scratch "$task.scratch"
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.7.4-foss-2016b'
  input:
    each cds from cds_path
    val nsample from adownsample

  output:
    path "${filename}" into log_counts
    
  script:
    filename = cds.getName()
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(leiden)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cluster_cds <- quickCluster(cds)
    cds <- computeSumFactors(cds, cluster=cluster_cds, min.mean=0.1)
    cds <- logNormCounts(cds, downsample=${nsample})
    saveRDS(cds, "${filename}")
    """  
}
// Get highly variable genes
process OSCA_SELECT{
  echo false
  scratch "$task.scratch"
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.7.4-foss-2016b'
  input:
    val cds from log_counts
    val nfeature from anfeatures

  output:
    path "${filename}" into osca_select
    
  script:
    filename = cds.getName()
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(leiden)
    set.seed(12357)

    cds <- readRDS("${cds}")
    gene_var <- modelGeneVar(cds)
    variable_genes <- getTopHVGs(gene_var, n=${nfeature})
    cds <- runPCA(cds, subset_row=variable_genes)
    pcs <- reducedDim(cds)
    choices <- getClusteredPCs(pcs)
    reducedDim(cds, "PCA_chosen") <- reducedDim(cds, "PCA")[,1:metadata(choices)\$chosen]
    metadata(cds)\$modelGeneVar <- gene_var
    saveRDS(cds, "${filename}")
    """  
}
// Duplicate input channels
osca_select.into{osca_combine; osca_data}
// Reduce dimensionality 
process OSCA_BATCHCORRECT {
  echo false
  scratch "$task.scratch"
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.7.4-foss-2016b'
  input:
    val sce_list from osca_combine.collect()
    val nfeature from anfeatures

  output:
    path "${study_name}_aggregated.rds" into merged_sce

  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(leiden)
    set.seed(12357)

    sce_files <- c("${sce_list.join('\",\"')}")
    # Get gene variation model
    getVar <- function(sce) {
      gene_var <- metadata(sce)\$modelGeneVar
      return(gene_var)
    }
    var_genes <- sce_files %>% map(getVar) 
    gene_comb <- do.call(combineVar, gene_var)
    gene_hvgs <- gene_comb$bio > 0
    sce_files <- do.call(multiBatchNorm, sce_files)
    merged <- do.call(cbind, sce_files)
    if (length(unique(merged\$Batch)) > 1) {
      merged <- fastMNN(merged, batch=merged\$Batch, 
                         subset.row = gene_hvgs, BSPARAM=BiocSingular::RandomParam(deferred = TRUE))
    }

    metadata(merged)\$Sample <- "${study_name}_aggregated"

    #Combine metdata
    combine_metrics <- function(sce) {
      sample <- metadata(sce)\$Sample
      metrics <- metadata(sce)\$perCellQCMetrics_filtered %>% as_tibble() %>% mutate(Sample = sample)
      return(metrics)
    }
    combine_info <- function(sce) {
      sample <- metadata(sce)\$Sample
      info <- metadata(sce)\$study_info
      return(info)
    }
    combine_vdj <- function(sce) {
      sample <- metadata(sce)\$Sample
      if (metadata(sce)\$vdj_raw != NULL) {
        metadata(sce)\$vdj_raw <- metadata(sce)\$vdj_raw %>% mutate(Sample = sample)
      }
      info <- metadata(sce)\$study_info
      return(info)
    }
    metadata(merged)\$perCellQCMetrics <- cds_list %>% map(combine_metrics) %>% purrr::reduce(rbind)
    metadata(merged)\$study_info <- cds_list %>% map(combine_info) %>% purrr::reduce(rbind)
    metadata(merged)\$vdj_raw <- cds_list %>% map(combine_vdj) %>% purrr::reduce(rbind)        
    saveRDS(merged, "${study_name}_aggregated.rds")
    """  
}
//Combine aggregated and sample channels
sce_path = osca_data.mix( merged_sce )
//Duplicate dimred channel
rmethod.into{ cluster_rmethod; plot_rmethod}
process OSCA_CLUSTER {
  echo false
  scratch "$task.scratch"
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.7.4-foss-2016b'
  input:
    each sce from sce_path
    val nsample from adownsample
    each nfeature from anfeatures
    each cnumber from acenters
    val method from cluster_rmethod
    val cluster from cmethod

  output:
    file "${filename}" into sce_cluster
    
  script:
    filename = sce.getName()
    sample = sce.getSimpleName()
    """ 
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(leiden)
    set.seed(12357)
    sce <- readRDS(${sce})
    if ("${method}"" == "TSNE"){
      sce <- runTSNE(sce, dimred="PCA_chosen") 
    } else {
      sce <- runUMAP(sce, dimred="PCA_chosen")
    }
    cluster <- buildSNNGraph(sce, k=${cnumber}, use.dimred = 'PCA')
    if ("${cluster}" == "leiden") {
      cluster_adj <- igraph::as_adjacency_matrix(cluster)
      cluster_leiden <- leiden(cluster_adj)
      sce\$cluster <- factor(cluster_leiden)
    } else {
      cluster_louvain <- igraph::cluster_louvain(cluster)
      sce\$cluster <- factor(cluster_louvain\$membership)
    }
    sce\$cluster <- factor(cluster_leiden)
    abi_run_condition <- tibble(Sample = metadata\$Sample, preprocess_num_dim = metadata(choices)\$chosen, 
                            preprocess_norm_method = "Log normal", preprocess_method = "PCA", 
                            reduction_method = "${method}", 
                            k = ${cnumber}, cluster_method = "${cluster}")
    metadata(sce) <- list(run_condition = abi_run_condition)
    saveRDS(sce, "${filename}_osca.rds")
    """
}
//Duplicate channels
sce_cluster.into{ sce_gather; sce_plot; sce_write}
//Plot UMAP or tSNE
process OSCA_ABINTIO_PLOT {
  echo false
  scratch "$task.scratch"
  publishDir "$params.output.folder/OSCA/Figures/" , mode : 'move'
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.7.4-foss-2016b'
  input:
    each sce from sce_cluster
    val rmethod from plot_rmethod

  output:
    file "${filename}_${rmethod}.png" into sce_plot
    
  script:
    filename = sce.getName()
    sample = sce.getSimpleName()
    """ 
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(leiden)
    set.seed(12357)

    sce <- readRDS("${sce}")
    plot <- plotReducedDim(sce, dimred="${rmethod}", colour_by="cluster")
    ggsave("${filename}_${rmethod}.png", plot, units = "in", height = 10, width = 10, dpi = "retina")
    """
}

process OSCA_GATHER {

  echo false
  publishDir "$params.output.folder/OSCA/Cluster/Metadata" , mode : 'move'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val sce_list from abinitio_sce.collect()
    val sid from study_id
  output:
    path "${sid}_metadata.csv" into sce_gather
  
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(DropletUtils)
  set.seed(12357)
  sce_list <- c("${sce_list.join('\",\"')}")
  getClusterSummary <- function(sce_file) {
    sce <- readRDS(sce_file)
    cluster_summary <- metadata(sce)\$run_condition
    return(cluster_summary)
  }
  cluster_summary <- sce_list %>% map(getClusterSummary) %>% reduce(rbind)
  write_csv(cluster_summary, "${sid}_metadata.csv")
  """
}

process OSCA_WRITE {

  echo false
  publishDir "$params.output.folder/OSCA/Cluster/SCE" , mode : 'move'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val sce from sce_write

  output:
    path "${file_name}" into sce_output
  
  script:
  file_name = sce.getName()
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(DropletUtils)
  set.seed(12357)

  sce <- readRDS("${sce}")
  writeRDS(sce, "${file_name}")
  """
}


/*process OSCA_PRIORI_CLUSTER {
  echo false
  publishDir "$params.output.folder/OSCA/Cluster/Priori"
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.7.4-foss-2016b'
  input:
    each cds from pr_cds_path
    val nsample from pdownsample
    each cnumber from pcenters

  output:
    file "Priori_${uuid}.rds" into priori_cds
    file "Priori_${uuid}_tSNE.png" into priori_tsne
    file "Priori_${uuid}_UMAP.png" into priori_umap
  
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(msigdbr)
    library(leiden)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cluster_cds <- quickCluster(cds)
    cds <- computeSumFactors(cds, cluster=cluster_cds, min.mean=0.1)
    cds <- logNormCounts(cds, downsample=${nsample})
    gene_var <- modelGeneVar(cds)
    #gene_var_cv2 <- modelGeneCV2(cds)
    immuno_sig <- msigdbr(species = "Homo sapiens", category = "C7")
    immuno_sig_genes <- rowData(cds)\$Symbol %in% immuno_sig\$human_gene_symbol
    cds <- runPCA(cds, subset_row=immuno_sig_genes)
    pcs <- reducedDim(cds)
    choices <- getClusteredPCs(pcs)
    reducedDim(cds, "PCA_chosen") <- reducedDim(cds, "PCA")[,1:metadata(choices)\$chosen]
    cds <- runTSNE(cds, dimred="PCA_chosen") 
    cds <- runUMAP(cds, dimred="PCA_chosen")
    cluster <- buildSNNGraph(cds, k=${cnumber}, use.dimred = 'PCA')
    cluster_adj <- igraph::as_adjacency_matrix(cluster)
    cluster_leiden <- leiden(cluster_adj)
    cluster_louvain <- igraph::cluster_louvain(cluster)
    #cds\$cluster <- factor(cluster_louvain\$membership)
    cds\$cluster <- factor(cluster_leiden)
    pri_run_condition <- tibble(uuid = "${uuid}", cds_file = paste(paste("Priori", "${uuid}", sep="_"), "rds", sep="."), preprocess_num_dim = metadata(choices)\$chosen, 
                            preprocess_norm_method = "Log normal", preprocess_method = "PCA", 
                            reduction_method = "UMAP/tSNE", 
                            k = ${cnumber}, cluster_method = "louvain", sample = "${sample}")
    metadata(cds) <- list(run_condition = pri_run_condition)
    saveRDS(cds, paste(paste("Priori", "${uuid}", sep="_"), "rds", sep="."))
    priori_tsne <- plotReducedDim(cds, dimred="TSNE", colour_by="cluster")
    priori_umap <- plotReducedDim(cds, dimred="UMAP", colour_by="cluster")
    ggsave(paste(paste("Priori", "${uuid}", "tSNE", sep="_"), "png", sep="."), priori_tsne, units = "in", dpi = "retina")
    ggsave(paste(paste("Priori", "${uuid}", "UMAP", sep="_"), "png", sep="."), priori_umap, units = "in", dpi = "retina")
    """

}

process PRI_GATHER_METADATA {

  echo false
  publishDir "$params.output.folder/OSCA/Cluster/Priori" , mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val cds_list from priori_cds.collect()
  output:
    path "Priori_cluster_metadata.csv" into pri_report
  
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(DropletUtils)
  set.seed(12357)
  cds_list <- c("${cds_list.join('\",\"')}")
  getClusterSummary <- function(cds_file) {
    cds <- readRDS(cds_file)
    cluster_summary <- metadata(cds)\$run_condition
    return(cluster_summary)
  }
  cluster_summary <- cds_list %>% map(getClusterSummary) %>% reduce(rbind)
  write_csv(cluster_summary, "Priori_cluster_metadata.csv")
  """
}


/*process GAT_UMAP_CLUSTER {
  echo false
  publishDir "$params.output.folder/Gattardo/Cluster/UMAPOpt"
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    each cds from um_cds_path
    val nsample from udownsample


  output:
    file "Monocle_${uuid}.rds" into monocle_cds
    
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      

    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(msigdbr)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cluster_cds <- quickCluster(cds)
    cds <- computeSumFactors(cds, cluster=cluster_cds, min.mean=0.1)
    cds <- logNormCounts(cds, downsample=${nsample})
    gene_var <- modelGeneVar(cds)
    #gene_var_cv2 <- modelGeneCV2(cds)
    immuno_sig <- msigdbr(species = "Homo sapiens", category = "C7")
    immuno_sig_genes <- rowData(cds)\$Symbol %in% immuno_sig\$gs_name\$human_gene_symbol
    cds <- runPCA(cds, subset_row=immuno_sig_genes)
    pcs <- reducedDim(cds)
    choices <- getClusteredPCs(cds)
    reducedDim(cds, "PCA_chosen") <- reducedDim(cds, "PCA")[,1:metdata(choices)\$chosen]
    cds <- runUMAP(cds, dimred="PCA_chosen", n_neighbors = ${uneighbor}, min_dist = ${udist})
    cluster <- buildSNNGraph(cds, k=${cnumber}, use.dimred = 'PCA')
    cluster_louvain <- igraph::cluster_louvain(cluster)
    cds\$cluster <- factor(cluster_louvain)
    run_condition <- tibble(uuid = "${uuid}", cds_file = paste(paste("Priori", "${uuid}", sep="_"), "rds", sep="."), preprocess_num_dim = metdata(choices)\$chosen, 
                            preprocess_norm_method = "Log normal", preprocess_method = "PCA", 
                            reduction_method = "UMAP/tSNE", 
                            k = ${cnumber}, cluster_method = "louvain", sample = cell_metadata\$Sample[1])
    metadata(cds) <- list(run_condition = run_condition)
    saveRDS(cds, paste(paste("Priori", "${uuid}", sep="_"), "rds", sep="."))
    priori_tsne <- plotReducedDim(cds, dimred="TSNE", colour_by="cluster")
    priori_umap <- plotReducedDim(cds, dimred="UMAP", colour_by="cluster")
    ggsave(paste(paste("Priori", "${uuid}", "tSNE", sep="_"), "png", sep="."), priori_tsne, units = "in", dpi = "retina")
    ggsave(paste(paste("Priori", "${uuid}", "UMAP", sep="_"), "png", sep="."), priori_umap, units = "in", dpi = "retina")
    """

}

process GAT_TSNE_CLUSTER {
  echo false
  publishDir "$params.output.folder/Gattardo/Cluster/TSNEOpt"
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    each cds from ts_cds_path
    val core from num_cores
    each pdim from num_dim
    val pnorm from norm_method
    val preduce from dimred_method
    each rreduce from reduce_method
    each cmethod from clustering_method
    each cnumber from cluster_number

  output:
    file "Monocle_${uuid}.rds" into monocle_cds
    
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      

    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(msigdbr)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cluster_cds <- quickCluster(cds)
    cds <- computeSumFactors(cds, cluster=cluster_cds, min.mean=0.1)
    cds <- logNormCounts(cds, downsample=${nsample})
    gene_var <- modelGeneVar(cds)
    #gene_var_cv2 <- modelGeneCV2(cds)
    immuno_sig <- msigdbr(species = "Homo sapiens", category = "C7")
    immuno_sig_genes <- rowData(cds)\$Symbol %in% immuno_sig\$gs_name\$human_gene_symbol
    cds <- runPCA(cds, subset_row=immuno_sig_genes)
    pcs <- reducedDim(cds)
    choices <- getClusteredPCs(cds)
    reducedDim(cds, "PCA_chosen") <- reducedDim(cds, "PCA")[,1:metdata(choices)\$chosen]
    cds <- runTSNE(cds, dimred="PCA_chosen", perplexity = ${perplexity}) 
    """

}*/
