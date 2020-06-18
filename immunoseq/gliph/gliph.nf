trb_files = Channel.fromPath(params.input.trb_path)

process gliph_prepare {
  echo false
  publishDir "$params.output.folder/gliph/input/", mode : "copy"
  label 'gizmo'
  scratch "$task.scratch"
  module 'R/3.6.2-foss-2016b-fh1'
  
  input:
    path tcr_path from trb_files
  output:
    path "*.tsv" into trb_byLength
  script:
    """
    #!/usr/bin/env Rscript
    .libPaths(c("/home/sravisha/R/x86_64-pc-linux-gnu-library/3.6", .libPaths()))  
    library(tidyverse)
    library(devtools)
    library(randomcoloR)
    library(patchwork)
    library(ggalluvial)
    devtools::load_all("/home/sravisha/projects/LymphoSeq2")

    createGliphByLength <- function(sample_table) {
        sample <- sample_table %>% 
                  pull(repertoire_id) %>% 
                  unique()
        sample_table <- sample_table %>% 
                        mutate(length = nchar(junction_aa)) %>%
                        filter(length >= 7) %>%
                        group_by(length) %>%
                        group_split(keep=TRUE) %>%
                        map(~ writeGliph(.x, sample)) 

    }

    writeGliph <- function(sample_table, sample) {
        length <- sample_table\$length[1]
        out_path <- paste(paste(sample, length, sep="_"), "tsv", sep=".")
        sample_table <- sample_table %>% 
                        mutate(Sample = repertoire_id) %>% 
                        select(junction_aa, Sample) %>% 
                        rename(CDR3b = junction_aa)
        write_tsv(sample_table, out_path)
    }

    data_path <- c(\"${tcr_path}\")
    stable <- LymphoSeq2::readImmunoSeq(data_path)
    atable <- LymphoSeq2::productiveSeq(stable)
    atable <- atable %>% 
              group_by(repertoire_id) %>% 
              group_split(keep=TRUE) %>% 
              map(createGliphByLength)
    """
}


process gliph_analyze {
  echo false 
  publishDir "$params.output.folder/gliph/output/", mode: 'copy'
  label 'gizmo'
  scratch "$task.scratch"
  stageInMode 'copy'
  module "Perl"

  input:
    each tcr_length from trb_byLength.collect().flatten()
    
  output:
    path "${tcr_length.getName()}-*.txt" into gliph_clonotypes

  script:
    sample = tcr_length.getSimpleName()    
    """
    perl /home/sravisha/software/gliph/gliph/bin/gliph-group-discovery.pl --tcr ${tcr_length} 
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

    sendMail(to: 'sravisha@fredhutch.org', subject: 'GLIPH nextflow execution', body: msg)
}