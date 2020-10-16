#!/usr/bin/env nextflow

//params.reads = "$params.input.fastq_path/*/*_L00{1,2}_R{1,2}_00{1,2,3,4}.fastq.gz"
//fastq_path = Channel.fromSRA('ERP105213')
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
chrom_path = Channel.fromPath(params.input.chrom_path)
cutsite_path = Channel.fromPath(params.input.cutsite_path)
config_path = Channel.fromPath(params.input.config_path)

process unzipFastq {
    echo false
    publishDir "$params.output.folder/preFastqQC/${sample}", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "FastQC"
    input:
        each accession from file("$params.input.acc_list").readLines()
    
    output:
        tuple val(accession), path("${accession}.fastq") into unzip_out

    script:
        """
        fasterq-dump -e $task.cpus ${accession} 
        """
}

unzip_out.into{preqc_path; trim_path}

process preFastQC {
    echo false
    publishDir "$params.output.folder/preFastqQC/${sample}", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "FastQC"
    input:
        set val(sample), path(read_one) from preqc_path
    
    output:
        tuple val(sample), path("*") into preqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} 
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
        set val(sample), path(read_one) from trim_path
    
    output:
        tuple val(sample), path("${sample}_trimmed.fastq"), path("${sample}_*.txt") into trim_out

    script:
        """
        bbduk.sh in=${read_one} out=${sample}_trimmed.fastq \\
        stats=${sample}_stats.txt statscolumns=5 \\
        bhist=${sample}_bhist.txt qhist=${sample}_qhist.txt qchist=${sample}_qchist.txt \\
        aqhist=${sample}_aqhist.txt bqhist=${sample}_bqhist.txt lhist=${sample}_lhist.txt \\
        phist=${sample}_phist.txt gchist=${sample}_gchist.txt k=13 t=$task.cpus  \\
        qtrim=rl tossbrokenreads=t -Xmx8g
        """
}

trim_out.into{align_path; postqc_path}

process postFastQC {
    echo false
    publishDir "$params.output.folder/postFastqQC/${sample}", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "FastQC"
    input:
        set val(sample), path(read_one), path(stats_path) from postqc_path
    
    output:
        tuple val(sample), path("*") into postqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} 
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
        set val(sample), path(read_one), path(stats_path) from align_path
    
    output:
        tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") into align_out

    script:
        """
        bowtie2 --very-sensitive --dovetail --met-file ${sample}_bmet.txt -p $task.cpus \\
        -x $params.input.index_path -U ${read_one} -S ${sample}.sam 
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
        tuple val(sample), path("${sample}.bam") into config_out

    script:
        """
        gatk SortSam -I=${sam_path} -O=${sample}_SR.bam -SO=coordinate --CREATE_INDEX=true -R=$params.input.fasta_path
        gatk MarkDuplicates -I=${sample}_SR.bam -M=${sample}_dmetrics.txt -O=${sample}.bam --CREATE_INDEX=true -R=$params.input.fasta_path
        gatk CollectAlignmentSummaryMetrics -I=${sample}.bam -O=${sample}_ametrics.txt
        """
    
}

process callPeaks { 
    echo false
    label 'gizmo'
    publishDir "$params.output.folder/peakCalls", mode : "copy"
    scratch "$task.scratch"
    module "MACS2"

    input:
        set val(sample), path(bam_path) from config_out
    output:
        tuple val("${sample}"), path("${sample}_treat_pileup.bdg"), path("${sample}_control_lambda.bdg") into bdg_out

    script:
        """
        macs2 callpeak -t ${bam_path} -B --SPMR --broad -g hs -n ${sample}
        """
}


process callFR {
    echo false
    publishDir "$params.output.folder/bedCompare", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "MACS2"

    input:
        set val(sample), path(treat), path(lambda) from bdg_out
    
    output:
        tuple val(sample), path("${sample}_FE.bdg"), path("${sample}_logLR.bdg") into bdgc_out

    script:
        """
        macs2 bdgcmp -t ${treat} -c ${lambda} -o ${sample}_FE.bdg -m FE
        macs2 bdgcmp -t ${treat} -c ${lambda} -o ${sample}_logLR.bdg -m logLR -p 0.00001
        """
}

process creatBigWig {
    echo false
    publishDir "$params.output.folder/bigWigs", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "bedtools"
    module "ucsc"

    input:
        set val(sample), path(fe_path), path(llr_path) from bdgc_out
        each cl_path from chrom_path
    
    output:
        tuple val(sample), path("${sample}_FE.bw"), path("${sample}_logLR.bw") into bwig_out

    script:
        """
        bedtools slop -i ${fe_path} -g ${cl_path} -b 0 > ${fe_path}.genome 
        bedClip ${fe_path}.genome ${fe_path} ${fe_path}.clip
        LC_COLLATE=C sort -k1,1 -k2,2n ${fe_path}.clip > ${fe_path}.sort.clip
        bedGraphToBigWig ${fe_path}.sort.clip ${cl_path} ${sample}_FE.bw
        rm -f ${fe_path}.clip ${fe_path}.sort.clip

        bedtools slop -i ${llr_path} -g ${cl_path} -b 0 > ${llr_path}.genome 
        bedClip ${llr_path}.genome ${llr_path} ${llr_path}.clip
        LC_COLLATE=C sort -k1,1 -k2,2n ${llr_path}.clip > ${llr_path}.sort.clip
        bedGraphToBigWig ${llr_path}.sort.clip ${cl_path} ${sample}_logLR.bw
        rm -f ${llr_path}.clip ${llr_path}.sort.clip
        """
}



process filterBigWig {
    echo false
    publishDir "$params.output.folder/finalBW", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "R/3.6.2-foss-2016b-fh1"

    input:
        set val(sample), path(fe_path), path(llr_path) from bwig_out
        each cut_site from cutsite_path
    
    output:
        tuple val(sample), path("${fe_path.getSimpleName()}_cutSite.bw"), path("${llr_path.getSimpleName()}_cutSite.bw") into fwig_out

    script:
        """
        #!/usr/bin/env Rscript
        library(rtracklayer)
        bigWigs <- import(\"${fe_path}\", format="BigWig")
        bedFile <- import(\"${cut_site}\", format="BED")
        seqlengths(bedFile) <- seqlengths(bigWigs)[names(seqlengths(bedFile))]
        overlap <- subsetByOverlaps(bedFile, bigWigs, type="any") 
        export(overlap, \"${fe_path.getSimpleName()}_cutSite.bw\", format="BigWig")
        bigWigs <- import(\"${llr_path}\", format="BigWig")
        bedFile <- import(\"${cut_site}\", format="BED")
        seqlengths(bedFile) <- seqlengths(bigWigs)[names(seqlengths(bedFile))]
        overlap <- subsetByOverlaps(bedFile, bigWigs, type="any") 
        export(overlap, \"${llr_path.getSimpleName()}_cutSite.bw\", format="BigWig")
        """
}



/*workflow.onComplete {

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

    sendMail(to: 'sravisha@fredhutch.org', subject: 'CutNRun nextflow execution', body: msg)
}*/
