trb_files = Channel.fromPath(params.input.trb_path + "/*_CFAR.tsv")
meta_file = Channel.fromPath(params.input.meta_path)

process network_prepare {
  echo false
  publishDir "$params.output.output_folder", mode : "copy"
  label 'gizmo_largenode'
  module 'R/4.0.2-foss-2019b'
  
  input:
    each tcr_path from trb_files
    val meta_data from meta_file
  output:
    path "${sample}_edist.tsv" into edit_path
  script:
    sample = tcr_path.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(LymphoSeq2)
    library(readxl)
    library(tidyverse)
    
    stable <- LymphoSeq2::readImmunoSeq(\"${tcr_path}\") %>%
              tidyr::uncount(duplicate_count) %>% 
              dplyr::mutate(duplicate_count = 1) 
    if (nrow(stable) < 10000) {
      sample <- nrow(stable)
    } else {
      sample <- 10000
    }
    stable <- stable %>%
              dplyr::sample_n(sample)
    atable <- LymphoSeq2::productiveSeq(stable) 
    ntable <- tidyr::expand_grid(amino_one = atable\$junction_aa, amino_two = atable\$junction_aa) %>% 
              dplyr::filter_all(all_vars(!is.na(.))) %>%
              dplyr::rowwise() %>%
              dplyr::mutate(editDist = adist(amino_one, amino_two, fixed=TRUE, partial=FALSE)) %>%
              dplyr::filter(editDist >= 1 & editDist <= 2) 
    sample <- stable %>% pull(repertoire_id) %>% unique()
    write_tsv(ntable, \"${sample}_edist.tsv\")
    """
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: 'sravisha@fredhutch.org', subject: 'ImmunoSeq network analysis nextflow execution', body: msg)
}
