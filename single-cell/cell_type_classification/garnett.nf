#!/usr/bin/env nextflow

cds_path = Channel.fromPath(params.input.cds_path + '/*.rds')
nunknowns = Channel.from(params.garnett.nunknowns)
known_types = Channel.from(params.garnett.known_type)
classifier_path = Channel.from(params.garnett.classifier_path)
classifier_type = Channel.from(params.garnett.classifier_type)
//public_trb = Channel.fromPath(params.input.ks_public)
CLASS_CH = classifier_path.merge(classifier_type)

cds_path.into{cds_gar; cds_vdj}
process GAR_CLASSIFY {
  echo false
  publishDir "$params.output.folder/Monocle/Classify/CDS" , mode : 'copy'
  module 'R/4.0.2-foss-2019b'
  module 'monocle3/0.2.2-foss-2019b-R-4.0.2'
  label 'mid_mem'
  input:
    each cds from cds_gar
    each nunknown from nunknowns
    each known from known_types
    set val(classifier), val(ctype) from CLASS_CH
  output:
    path "${sample}_classifier.rds" into garnett_out
    path "${sample}_cds.rds" into garnett_cds
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    library(garnett)
    library(org.Hs.eg.db)
    set.seed(12357)

    cds <- readRDS("${cds}")
    gene_metadata <- rowData(cds)
    cell_metadata <- colData(cds)
    if (${known} == TRUE) {
      cds\$known_type = metadata(cds)\$vdj_table  
    }
    marker_check <- check_markers(cds, "${classifier}", 
                                  db = org.Hs.eg.db, 
                                  cds_gene_id_type = "ENSEMBL", 
                                  marker_file_gene_id_type = "${ctype}")
    classifier <- train_cell_classifier(cds = cds, 
                                        marker_file = "${classifier}", 
                                        db = org.Hs.eg.db, 
                                        cds_gene_id_type = "ENSEMBL", 
                                        marker_file_gene_id_type = "${ctype}", 
                                        num_unknown = ${nunknown}, 
                                        cores = ${task.cpus})
    cds <- classify_cells(cds, 
                          classifier, db = org.Hs.eg.db, 
                          cluster_extend = TRUE, 
                          cds_gene_id_type = "ENSEMBL")               
    run_condition <- metadata(cds)\$run_condition %>%  
                     as_tibble() %>%
                     dplyr::mutate(gar_uuid = "${sample}", 
                                   gar_cds_file = paste(paste("Garnett", "${sample}", "cds", sep="_"), "rds", sep="."), 
                                   classifier = "${classifier}", num_unknowns = ${nunknown})
    metadata(cds) <- list(run_condition = run_condition)
    saveRDS(cds, paste(paste("${sample}", "cds", sep="_"), "rds", sep="."))
    saveRDS(classifier, paste(paste("${sample}", "classifier", sep = "_"), "rds", sep = "."))
    """
}
/*
process VDJ_UMAPS{
  echo false
  publishDir "$params.output.folder/Monocle/Classify/Plot" , mode : 'copy', pattern : '*.png'
  publishDir "$params.output.folder/Monocle/Classify/Tables" , mode : 'copy', pattern : '*.tsv'
  scratch "$task.scratch"
  module 'R/4.0.2-foss-2019b'
  module 'monocle3/0.2.2-foss-2019b-R-4.0.2'
  label 'low_mem'
  input:
    each cds from cds_vdj
    path ptable from public_trb
  output:
    path "${sample}_publicTRB.png" into vdj_plot
    path "${sample}_clonotypes.tsv" into vdj_tables
  script:
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    library(garnett)
    library(tidyverse)
    library(org.Hs.eg.db)
    set.seed(12357)

    cds <- readRDS("${cds}")
    ptable <- read_csv("${ptable}") %>% pull(aminoAcid)
    barcodes <- rownames(reducedDims(cds)[["UMAP"]])
    reducedDims(cds)[["UMAP"]] %>% as_tibble() %>% print()
    coord_table <- reducedDims(cds)[["UMAP"]] %>%
                   as_tibble() %>%
                   dplyr::rename(UMAP1 = V1,
                          UMAP2 = V2) %>%
                   dplyr::mutate(barcode = barcodes,
                          id_barcode = str_extract(barcode, "[ACGT]+"),
                          id_number = str_extract(barcode, "\\\\d+"))
    vdj_table <- metadata(cds)\$vdj_table
    sample <- metadata(cds)\$Sample
    vdj_pal <- c("CASSVDKGGTDTQYF" = "#a6cee3", "CASSLGQGSTDTQYF" = "#1f78b4",
                 "CASSVVGGTDTQYF" = "#b2df8a", "CASSIAGHEQYF" = "#33a02c",
                 "CASSLWGGPSNEQFF" = "#fb9a99", "CASSQEFPGAAYNEQFF" = "#e31a1c",
                 "CASSQEFPGALYNEQFF" = "#fdbf6f", "CASSQEFPGASYNEQFF" = "#ff7f00",
                 "CASSYPQGMNTEAFF" = "#cab2d6", "Absent" = alpha("grey", 0.4))
    coord_table <- dplyr::left_join(coord_table, vdj_table, by = c("Sample" = "repertoire_id", "barcode" = "Barcode")) %>%
                   dplyr::select(barcode, UMAP1, UMAP2, cdr3, chain) %>%
                   tidyr::separate_rows(cdr3, chain, sep = ";") %>%
                   dplyr::filter(chain == "TRB") %>%
                   dplyr::mutate(public_seq = if_else(cdr3 %in% ptable, cdr3, "Absent"))
    vdj_plot <- ggplot(coord_table, aes(x=UMAP1, y=UMAP2, color = public_seq, shape = public_seq)) +
                geom_point() +
                scale_color_manual(values = vdj_pal) +
                theme_classic() +
                ggtitle(sample) +
                theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(vdj_plot, filename = "${sample}_publicTRB.png", height = 10, width = 12, units = "in", device = "png")
    write_tsv(vdj_table, "${sample}_clonotypes.tsv")
    """
}*/


garnett_cds.into{summary_cds; plot_cds}

process GAR_SUMMARIZE_RUN {
  echo false
  publishDir "$params.output.folder/Monocle/Classify/Metrics" , mode : 'copy'
  module 'R/4.0.2-foss-2019b'
  module 'monocle3/0.2.2-foss-2019b-R-4.0.2'
  label 'low_mem'
  input:
    path cds from summary_cds.collect()

  output:
    path "Run_summary.csv" into run_summary

  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    library(garnett)
    library(org.Hs.eg.db)
    set.seed(12357)

    cds_list <- list.files(path = ".", pattern = "*_cds.rds")

    summarizeCDS <- function(cds_file) {
      cds <- readRDS(cds_file)
      discard_summary <- metadata(cds)\$run_condition
      return(discard_summary)       
    }

    run_summary <- cds_list %>% map(summarizeCDS) %>% bind_rows()
    write_csv(run_summary, 'Run_summary.csv')
    """    
}

process GAR_SUMMARIZE_CLASSIFIER {
  echo false
  publishDir "$params.output.folder/Monocle/Classify/Metrics" , mode : 'copy'
  module 'R/4.0.2-foss-2019b'
  module 'monocle3/0.2.2-foss-2019b-R-4.0.2'
  label 'low_mem'
  input:
    path classifier from garnett_out.collect()

  output:
    file "Classifier_summary.csv" into classifier_summary

  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    library(garnett)
    library(org.Hs.eg.db)
    set.seed(12357)

    classifier <- list.files(path = ".", pattern = "*_classifier.rds")

    summarizeClassifier <- function(cds_file){
      cds <- readRDS(cds_file)
      uuids <- str_split(cds_file, "_")[1]
      cell_type_summary <- colData(cds) %>% 
                           as_tibble() %>% 
                           dplyr::mutate(uuid = uuids) %>% 
                           select(uuid, Sample, Barcode, cell_type)
      
    }
    classifier_summary <- classifier %>% 
                          map(summarizeClassifier) %>% 
                          bind_rows()
    classifier_summary <- classifier_summary %>% 
                          group_by(sample, uuid, cell_type) %>% 
                          summarize(count = n()) %>% 
                          pivot_wider(id_cols = sample, 
                                      names_from = cell_type, 
                                      values_from = count)

    write_csv(classifier_summary, 'Garnett_summary.csv')
    """    
}

process GAR_PLOT {
  echo false
  publishDir "$params.output.folder/Monocle/Classify/Plot" , mode : 'copy'
  module 'R/4.0.2-foss-2019b'
  module 'monocle3/0.2.2-foss-2019b-R-4.0.2'
  label 'mid_mem'
  input:
    each cds from plot_cds

  output:
    path "${sample}_cell_classification.png" into plot_out
    

  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    library(garnett)
    library(org.Hs.eg.db)
    set.seed(12357)

    cds <- readRDS("${cds}")
    sample <- "${sample}"
    icr_palate <- c("Naive CD4 T cells" = "#f54242", "T regulatory cells" = "#eb2121", "Follicular helper T cells" = "#c24242",
    "Th1 cells" = "#8c2323", "Th1/Th17 cells" = "#7a0909", "Th17 cells" = "#b81f1f" , "Th2 cells" = "#c90808",
    "Helper T cells" = "#e81010", "Follicular helper T cells and Th1/Th17 cells" = "#cf0e0e",
    "Naive CD4 T cells and T regulatory cells" = "#9e0606", "Naive T helper cells" = "#f5140c",
    "Follicular helper T cells and Helper T cells" = "#f54c4c", "Terminal effector CD4 T cells" = "#f20c0c", "Naive CD4 T cells" = "#f5140c",
    "Memory CD4 T cells" = "#fa0202", "CD4 T cells" = "#de0202", "Naive CD8 T cells" = "#3dfa1b", "Central memory CD8 T cells" = "#2fe60e",
    "Effector memory CD8 T cells" = "#2bc90e", "Terminal effector CD8 T cells" = "#24ad0a", "T cells" = "#30bf00",
    "Effector memory and Terminal effector CD8 T cells" = "#1e8f09", "Memory CD8 T cells" = "#168701", "CD8 T cells" = "#127500",
    "Naive B cells" = "#4f8ce8", "Non-switched memory B cells" = "#4476c2", "Switched memory B cells" = "#305da1",
    "Exhausted B cells" = "#2373eb", "Exhausted and Switched memory B cells" = "#0d5dd6",
    "Naive and Non-switched memory B cells" = "#0c4bab", "Naive B cells" = "#156bed", "Memory B cells" = "#18478f",
    "Naive and Memory B cells" = "#0448b0", "Plasmablasts" = "#023178", "B cells" = "#004fc7", "Plasmacytoid dendritic cells" = "#004fc7",
    "Myeloid dendritic cells" = "gold1", "Classical monocytes" = "gold2", "Intermediate monocytes" = "gold3",
    "Non-classical monocytes" = "#663605", "Non-classical and Intermediate monocytes" = "#6b3c0c", "Monocytes" = "#7d4002",
    "Dendritic cells" = "#734c27", "Myeloid Phagocytes" = "#804208", "Vd2 gamma delta T cells" = "#d145b3",
    "Non Vd2 gamma delta T cells" = "#db2ab5", "Gamma delta T cells" = "#ed3bc7", "MAIT cells" = "#d406a8", "Innate T cells" = "#f507c2",
    "Natural killer cells" = "#00f7f7", "Low density neutrophils" = "#ff8400", "Low density basophils" = "#f5850c",
    "Low density granulocytes" = "#fa8f1b", "Progenitor cells" = "#f5f502", "Unknown" = "#bfbfbd", "Myeloid" = "#4a2e01")
    dim_name <- reducedDimNames(cds)
    if ("UMAP" %in% dim_name ) {
      umap_plot <- plot_cells(cds, reduction_method="UMAP", group_cells_by="partition", color_cells_by="cluster_ext_type", label_cell_groups=FALSE, cell_size=0.5) + scale_color_manual(values=icr_palate)
      ggsave(paste(sample, "cell_classification.png", sep="_"), plot=umap_plot, height=14.11, width=28.22, units="cm", dpi="retina")
    }
    if ("tSNE" %in% dim_name) {
      tsne_plot <- plot_cells(cds, reduction_method="tSNE", group_cells_by="partition", color_cells_by="cluster_ext_type", label_cell_groups=FALSE, cell_size=0.5) + scale_color_manual(values=icr_palate)
      ggsave(paste(sample, "cell_classification.png", sep="_"), plot=tsne_plot, height=14.11, width=28.22, units="cm", dpi="retina")
    }
    """    
}
