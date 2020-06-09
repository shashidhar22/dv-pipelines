#!/usr/bin/env nextflow

params.reads = "$params.input.fastq_path/*/*_L00{1,2}_R{1,2}_00{1,2,3,4}.fastq.gz"
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
config_path = Channel.fromPath(params.input.config_path)

process combineFastq {
    echo false
    publishDir "$params.output.folder/combinedFastq", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    input:
        set sample, path(fastq_group) from fastq_path

    output:
        tuple val(sample), path("${sample}_R1.fastq"), path("${sample}_R2.fastq") into comb_out

    script:
        """
        zcat *R1*.fastq.gz > ${sample}_R1.fastq
        zcat *R2*.fastq.gz > ${sample}_R2.fastq
        """ 
    }

comb_out.into{preqc_path; trim_path}

process preFastQC {
    echo false
    publishDir "$params.output.folder/preFastqQC/${sample}", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "FastQC"
    input:
        set val(sample), path(read_one), path(read_two) from preqc_path
    
    output:
        tuple val(sample), path("*") into preqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
        """
}

process trimFastq {
    echo false
    publishDir "$params.output.folder/trimFastq/${sample}", pattern: "*.fastq", mode : "copy"
    publishDir "$params.output.folder/trimFastq/${sample}/Stats", pattern: "*.txt", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "BBMap"
    input:
        set val(sample), path(read_one), path(read_two) from trim_path
    
    output:
        tuple val(sample), path("${sample}_trimmed_R1.fastq"), path("${sample}_trimmed_R2.fastq"), path("${sample}_*.txt") into trim_out

    script:
        """
        bbduk.sh in=${read_one} in2=${read_two} out=${sample}_trimmed_R1.fastq \\
        out2=${sample}_trimmed_R2.fastq stats=${sample}_stats.txt statscolumns=5 \\
        bhist=${sample}_bhist.txt qhist=${sample}_qhist.txt qchist=${sample}_qchist.txt \\
        aqhist=${sample}_aqhist.txt bqhist=${sample}_bqhist.txt lhist=${sample}_lhist.txt \\
        phist=${sample}_phist.txt gchist=${sample}_gchist.txt k=13 t=$task.cpus  \\
        qtrim=rl -Xmx8g
        """
}

trim_out.into{align_path; postqc_path}

process postFastQC {
    echo false
    publishDir "$params.output.folder/postFastqQC/${sample}", mode : "move"
    label 'gizmo'
    scratch "$task.scratch"
    module "FastQC"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from postqc_path
    
    output:
        tuple val(sample), path("*") into postqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
        """
}

process alignReads {
    echo false
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.sam", mode : "copy"
    publishDir "$params.output.folder/alignedReads/${sample}/Stats", pattern: "*.txt", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "Bowtie2"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from align_path
    
    output:
        tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") into align_out

    script:
        """
        bowtie2 --very-sensitive --dovetail --met-file ${sample}_bmet.txt -p $task.cpus \\
        -x $params.input.index_path -1 ${read_one} -2 ${read_two} -S ${sample}.sam 
        """ 
}

process processAlignments {
    echo false
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.ba*", mode : "copy"
    publishDir "$params.output.folder/alignedReads/${sample}/Stats", pattern: "*.txt", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "GATK"
    input:
        set val(sample), path(sam_path), path(align_met) from align_out
    
    output:
        tuple val(sample), path("${sample}.bam"), path("${sample}.bai"), path("${sample}_*.txt") into postal_out
        path "${sample}.bam" into config_out
    script:
        group = ${sample} =~ /
        """
        gatk SortSam -I=${sam_path} -O=${sample}_SR.bam -SO=coordinate --CREATE_INDEX=true -R=$params.input.fasta_path
        gatk MarkDuplicates -I=${sample}_SR.bam -M=${sample}_dmetrics.txt -O=${sample}.bam --CREATE_INDEX=true -R=$params.input.fasta_path
        gatk CollectAlignmentSummaryMetrics -I=${sample}.bam -O=${sample}_ametrics.txt
        """
    
}

process callPeaks { //.splitCsv(header: true, quote: '\"')
    echo false
    publishDir "$params.output.folder/peakCalls/${sample}", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "MACS2"
    input:
        each sample from config_path.splitCsv(header: true, quote: '\"')
        val results from postal_out.collect()
    output:


    script:
        """
        macs2 callpeak -t $sample.condition -c $sample.control -B --SPMR --broad -g hs -n $sample.sample_$sample.cond_$sample.cont
        """
}

aling_out.group_by("$params.output.folder/alignedReads/*.bam", size:-1).println()
