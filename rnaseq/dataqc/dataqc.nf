fastq_path = Channel.fromFilePairs(params.input.fastq_path + '/*_R{1,2}*.fq.gz', size: 2)
bam_path = Channel.fromPath(params.input.bam_path + '/*/*.Aligned.bam')
ref_path = Channel.fromPath(params.input.ref_path)
gtf_path = Channel.fromPath(params.input.gtf_path)


fastq_path.into{ qc_path, align_path}
process FASTQC {
    echo false
    scratch "$task.scratch"
    publishDir "$params.output.folder/FastqQC" , mode : 'move'
    module "fastqc/0.11.3"
    label "gizmo_largenode"

    input:
    tuple sample, files from fastq_path
    
    output:
    path "*" into fastqc_out

    script:
    rone = files[0]
    rtwo = files[1]
    file_reg = rone =~ /\d+-\d+-\w-\d+/
    file_name = file_reg[0]
    """
    zcat ${rone} ${rtwo} > ${file_name}.fastq
    fastqc -t $task.cores ${file_name}.fastq
    """
}

process BOWTIE {
    echo false
    scratch "$task.scratch"
    publishDir "$params.output.folder/FastqQC" , mode : 'move'
    module "Bowtie2"
    label "gizmo_largenode"

    input:
    tuple sample, files from fastq_path
    
    output:
    path "*" into fastqc_out

    script:
    rone = files[0]
    rtwo = files[1]
    file_reg = rone =~ /\d+-\d+-\w-\d+/
    file_name = file_reg[0]
    """
    zcat ${rone} ${rtwo} > ${file_name}.fastq
    fastqc -t $task.cores ${file_name}.fastq
    """

}


/*process featureCounts {
    echo false
    scratch "$task.scratch"
    publishDir "$params.output.folder/featureCounts" , mode : 'move'
    module "R/3.6.1-foss-2016b-fh2"
    label "gizmo_largenode"

    input:
    val bam_file from bam_path.collect()
    val gtf_file from gtf_path
    
    output:
    path "KS_RNA_Seq_counts.csv" into counts_out

    script:
    """
    #!/usr/bin/env Rscript
    library(Rsubread)
    library(tidyverse)

    bam_list <- c("${bam_file.join('\",\"')}")
    fc <- featureCounts(bam_list, annot.ext = "${gtf_file}", isGTFAnnotationFile = TRUE)
    write_csv(fc, "KS_RNA_Seq_counts.csv")
    """
}*
