trb_files = Channel.fromPath(params.input.trb_path)
meta_file = Channel.fromPath(params.input.meta_path)

trb_files.into{gliph_files; net_files; deep_files; sonia_files}

process gliph_prepare {
  echo false
  publishDir "$params.output.interm_folder", mode : "copy"
  label 'local'
  module 'R/4.0.2-foss-2019b'
  
  input:
    val tcr_path from gliph_files
    val meta_data from meta_file
  output:
    path "*.tsv" into study_trb
  script:
    """
    #!/usr/bin/env Rscript
    library(LymphoSeq2)
    library(readxl)
    library(tidyverse)
    
    createGliph <- function(atable) {
      patient_name <- atable %>%
                      pull(patientID) %>%
                      unique()
      out_file <- paste(patient_name, "tsv", sep = ".")
      gtable <- atable %>%
                select(junction_aa, v_call, j_call, sampleID) %>%
                rename(Sample = sampleID,
                       TRBV = v_call,
                       TRBJ = j_call,
                       CdR3b = junction_aa)
      write_tsv(gtable, out_file)
    }
    study_table <- LymphoSeq2::readImmunoSeq(\"${tcr_path}\")
    meta_table <- readxl::read_excel(\"${meta_data}\", sheet = "CFAR_Sample_info")
    meta_table <- meta_table %>%
                  group_by(patientID) %>%
                  arrange(dateTCR) %>%
                  filter(!is.na(dateTCR)) %>%
                  filter(patientID == "1218")      
    atable <- LymphoSeq2::productiveSeq(study_table) 
    atable <- left_join(meta_table, atable, by = c("sampleID" = "repertoire_id"))
    atable %>% 
    group_by(patientID) %>%
    group_split() %>%
    map(createGliph)
    """
}


process gliph_analyze {
  echo false 
  publishDir "$params.output.output_folder", mode: 'copy'
  label 'gizmo'
  scratch "$task.scratch"
  stageInMode 'copy'
  module "Perl"

  input:
    path patient from study_trb
    
  output:
    tuple val(sample), path("${patient.getName()}-convergence-groups.txt") into gliph_clonotypes

  script:
    sample = patient.getSimpleName() 
    """
    perl /home/sravisha/software/gliph/gliph/bin/gliph-group-discovery.pl --tcr ${patient} 
    """

}

/*process gliphCombine {
  echo false 
  publishDir "$params.output.folder/gliphByPT/combine/${sample}", mode: 'copy'
  label 'gizmo'
  scratch "$task.scratch"
  stageInMode 'copy'
  module "Perl"

  input:
    set val(sample), path(tcr_paths) from gliph_clonotypes.groupTuple()
  output:
    path "${sample}-convergence-groups.txt" into gliph_combined

  script:
    """
    cat *-convergence-groups.txt > ${sample}-convergence-groups.txt
    """
}*/


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
