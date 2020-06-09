#!/usr/bin/env nextflow

fastq_path = Channel.fromFilePairs('/home/sravisha/fngs/tenX/MMRFINE_10Xv3/Fastq/VDJ/*_L{001,002}_R{1,2}*.fastq.gz', size: 4)
ref_path = Channel.fromPath(params.input.ref_path)

process TRIM_READS {
  echo true
  //scratch "/fh/scratch/delete30/warren_h/sravisha/"
  //publishDir "$params.output.folder/ImmuneInsertions/TrimmedReads" , mode : 'copy'
  module 'BBMap'
  input:
    each fastq from fastq_path
    path ref from ref_path

  //output:
   // file "${sample}_trimmed_R1.fastq" into rones
   // file "${sample}_trimmed_R2.fastq" into rtwos

  script:

    """
    echo ${fastq}
    """
}

/*process BWA_ALING {
  echo false
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  publishDir "$params.output.folder/ImmuneInsertions/Alignment" , mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'
  label 'gizmo'
  input:
    each tfastq from trimmed_fastq
    path ref from ref_path

  output:
    file "${sample}_aligned.bam" into aligned

  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      
    sample = cds.getSimpleName()
    """
    """
}*/

