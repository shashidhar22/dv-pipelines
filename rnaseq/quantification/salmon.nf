#!/usr/bin/env nextflow

params.reads = "$params.input.fastq_path/*R{1,2}.fq.gz"
fastq_path = Channel.fromFilePairs(params.reads)
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
gff_path = Channel.fromPath(params.input.gff_path)


//Index references
process salmonIndex {
  echo false
  publishDir "$params.output.index", mode : "copy"
  label 'gizmo'
  scratch "$task.scratch"
  input:
    path ref_file from fasta_path
    
  output:
    path "${ref_file.getSimpleName()}" into index_paths

  script:
    """
    ~/software/salmon-latest_linux_x86_64/bin/salmon index -t ${ref_file} -i ${ref_file.getSimpleName()}
    """ 
 
}

// Run salmon
process salmonQuant {
  echo false
  publishDir "$params.output.folder", mode : "copy"
  label 'gizmo'
  scratch "$task.scratch"
  input:
    each fastq_data from fastq_path
    path index_file from index_paths 
    
  output:
    path "${sample}_quant" into meta_obj

  script:
    sample = fastq_data[0]
    fastq_files = fastq_data[1]
    """
    ~/software/salmon-latest_linux_x86_64/bin/salmon quant -i ${index_file}  -l A -1 ${fastq_files[0]} -2 ${fastq_files[1]} -p 4 --validateMappings -o ${sample}_quant
    """
    

}


//Generate count matrix 
/*process countMatrix {
  echo false
  publishDir "$params.output.folder", mode : "copy"
  label 'gizmo'
  scratch "$task.scratch"
  module "R/3.6.2-foss-2016b"
  input:
    val results_files from meta_obj.collect()
    val gff from gff_path
    val outdir from out_path
    
  output:
    path "KS_RNA_Gene_level.csv" into gene_matrix

  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(org.Hs.eg.db)
    library(tximport)
    library(biomaRt)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(ensembldb)
    library(EnsDb.Hsapiens.v86)

    ##Gencode analysis
    folder_list <- list.files(path = "/home/sravisha/fuser/kstme/RNASeq/quants", 
                              full.names = TRUE, 
                              recursive = TRUE, 
                              pattern = "quant.sf")
    meta_table <- tibble(samples = dirname(folder_list), 
                        paths = folder_list, 
                        genome = "Ensembl HG28", 
                        analysis = "Salmon")
    meta_table <- meta_table %>% 
                  mutate(samples = str_extract(samples, "008-\\d+-\\w-\\d+")) %>% 
                  mutate(samples = str_replace_all(samples, "-", "_")) %>%
                  mutate(samples = paste("S", samples, sep=""))
    sample_list <- meta_table %>% 
                  pull(paths)
    names(sample_list) <- meta_table %>% 
                          pull(samples)
    txdb <- makeTxDbFromGFF("/home/sravisha/fngs/ReferenceGenomes/Human_genomes/Gencode/gencode.v34.annotation.gff3.gz", 
                            format="gff3", 
                            dataSource = "Gencode34", 
                            organism = "Homo sapiens" )
    key_list <- keys(txdb, 
                    keytype = "TXNAME")
    tx2gene <- select(txdb, key_list, "GENEID", "TXNAME")
    txi <- tximport(sample_list, 
                    type = "salmon", 
                    tx2gene = tx2gene)
    ensembl <-  useMart("ensembl",
                        dataset="hsapiens_gene_ensembl")
    study_table <- txi\$counts %>% 
                  as_tibble() %>% 
                  mutate(GeneSymbol = rownames(txi$\$counts)) %>% 
                  dplyr::select(GeneSymbol, everything())
    gene_list <- getBM(attributes=c('ensembl_gene_id_version', 'hgnc_symbol'), 
                      filters = 'ensembl_gene_id_version', 
                      values = study_table\$GeneSymbol, 
                      mart = ensembl, 
                      uniqueRows = FALSE) %>% 
                as_tibble()
    study_table <- inner_join(study_table, gene_list, 
                              by=c("GeneSymbol" = "ensembl_gene_id_version")) %>% 
                  dplyr::select(-GeneSymbol) %>% 
                  rename("hgnc_symbol" = "GeneSymbol") %>% 
                  dplyr::select(GeneSymbol, everything())
    write_csv(study_table, path = "KS_RNA_Gene_level.csv")
    """
}*/