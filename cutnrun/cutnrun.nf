#!/usr/bin/env nextflow

params.reads = "$params.input.fastq_path/*/*_L00{1,2}_R{1,2}_00{1,2,3,4}.fastq.gz"
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
chrom_path = Channel.fromPath(params.input.chrom_path)
cutsite_path = Channel.fromPath(params.input.cutsite_path)
config_path = Channel.fromPath(params.input.config_path)
spikein_path = Channel.fromPath(params.input.spikeins)
spikein_geno = Channel.from(params.input.spike_gen)
spikein_fast = Channel.fromPath(params.input.spike_fasta)
plot_group = Channel.from(params.input.plot_groups)

process combineFastq {
    echo false
    publishDir "$params.output.folder/combinedFastq", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    stageInMode "copy"
    module "SRA-Toolkit"
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

trim_out.into{spike_path; postqc_path}

process postFastQC {
    echo false
    publishDir "$params.output.folder/postFastqQC/${sample}", mode : "copy"
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

process spikeInAlignment {
    echo false
    publishDir "$params.output.folder/spikedAlign/alignemnts/${sample}", pattern: "*.sam", mode : "copy"
    publishDir "$params.output.folder/spikedAlign/alignments/${sample}/Stats", pattern: "*.txt", mode : "copy"
    publishDir "$params.output.folder/spikedAlign/sampleFastq/", pattern: "*.fastq", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "Bowtie2"
    module "GATK"
    module "SAMtools"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from spike_path

    output:
        tuple val(sample), path("${sample}.bam.bai"), path("${sample}.bam"), path("${sample}*.fastq") into spike_out
        path "${sample}*.fastq" into human_fastq
        path "${sample}.ba*" into bam_path
    script:
        """
        bowtie2 --very-sensitive -p $task.cpus --un-conc ${sample} \\
        -x $params.input.spikeins -1 ${read_one} -2 ${read_two} -S ${sample}.sam 
        samtools view -b -S -h -o ${sample}_in.bam ${sample}.sam
        samtools sort ${sample}_in.bam -o ${sample}.bam
        samtools index ${sample}.bam 
        mv ${sample}.1 ${sample}_R1.fastq
        mv ${sample}.2 ${sample}_R2.fastq
        """ 
}

spike_out.into{ calc_spike; align_path}
human_fastq.into{ human_calc; human_subset}


process calcSpikeAbundance {
    echo false
    label 'gizmo'
    scratch "$task.scratch"
    publishDir "$params.output.folder/spikedAlign/normalization", pattern: "*.csv", mode : "copy"
    module "Python/3.7.4-foss-2019b-fh1"
    stageInMode "copy"

    input:
        path bam_files from bam_path.collect()
        path fastq_files from human_calc.collect()
    output:
        path "CNR_config.csv" into spike_config
    script:
        $/
        #!/usr/bin/env python

        import pysam 
        import glob
        import numpy as np

        def getReadCount(fastq_file):
            count = 0
            with open(fastq_file) as read:
                for id in read:
                    count += 1
                    next(read)
                    next(read)
                    next(read)
            return(count)

        def getMapped(bam_file):
            bam_reader = pysam.AlignmentFile(bam_file, "rb")
            mapped = bam_reader.get_index_statistics()[0].mapped
            return(mapped)

        bam_files = glob.glob("*.bam")
        bam_mapped = np.array([getMapped(bfile) for bfile in bam_files])
        print(bam_mapped)
        bam_ratio = min(bam_mapped)/bam_mapped
        rone_path = glob.glob("*R1.fastq")
        rtwo_path = glob.glob("*R2.fastq")

        out_file = "CNR_config.csv"
        out_writer = open(out_file, 'w')
        out_writer.write("sample,rone,rtwo,rcount\n")
        count = 0
        for rone, rtwo in zip(rone_path, rtwo_path):
            sample = rone.split('_R1')[0]
            rone_count = getReadCount(rone)
            normalized = bam_ratio[count] * rone_count
            out_writer.write("{0},{1},{2},{3}\n".format(sample, rone, rtwo, int(normalized)))
            count += 1
        
        out_writer.close()
        /$

}

process spikeNormalize {
    echo false
    label 'gizmo'
    scratch "$task.scratch"
    publishDir "$params.output.folder/normalizedFastq", pattern : "*.fq", mode : "copy"
    module "seqtk"
    stageInMode "copy"
    input:
        each groups from spike_config.splitCsv(header: true, quote: '\"')
        path  fastq_files from human_subset.collect()
    output:
        tuple val("${groups.sample}"), path("*_normed_*.fq") into normed_out
    script:
        """
        seqtk sample -s100 $groups.rone $groups.rcount > ${groups.sample}_normed_R1.fq
        seqtk sample -s100 $groups.rtwo $groups.rcount > ${groups.sample}_normed_R2.fq
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
        set val(sample), path(fastq_files) from normed_out
    
    output:
        tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") into align_out

    script:
        read_one = fastq_files[0]
        read_two = fastq_files[1]
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
        tuple val(sample), path("${sample}.bam"), path("${sample}.bai") into plot_out

    script:
        """
        gatk SortSam -I ${sam_path} -O ${sample}_SR.bam -SO coordinate --CREATE_INDEX true -R $params.input.fasta_path
        gatk MarkDuplicates -I ${sample}_SR.bam -M ${sample}_dmetrics.txt -O ${sample}.bam --CREATE_INDEX true -R $params.input.fasta_path
        gatk CollectAlignmentSummaryMetrics -I ${sample}.bam -O ${sample}_ametrics.txt
        """
    
}

process callPeaks { 
    echo false
    label 'gizmo'
    publishDir "$params.output.folder/peakCalls", mode : "copy"
    scratch "$task.scratch"
    module "MACS2"

    input:
        each sample from config_path.splitCsv(header: true, quote: '\"')
        val results from postal_out.collect()
    output:
        tuple val("$sample.sample"+"_"+"$sample.cond"+"_"+"$sample.cont"), path("$sample.sample"+"_"+"$sample.cond"+"_"+"$sample.cont"+"_treat_pileup.bdg"), path("$sample.sample"+"_"+"$sample.cond"+"_"+"$sample.cont"+"_control_lambda.bdg") into bdg_out

    script:
        """
        macs2 callpeak -t $sample.condition -c $sample.control -B --SPMR --broad -g hs -n $sample.sample"_"$sample.cond"_"$sample.cont
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
    module "BEDTools"

    input:
        set val(sample), path(fe_path), path(llr_path) from bdgc_out
        each cl_path from chrom_path
    
    output:
        tuple val(sample), path("${sample}_FE.bw"), path("${sample}_logLR.bw") into bwig_out

    script:
        """
        bedtools slop -i ${fe_path} -g ${cl_path} -b 0 > ${fe_path}.genome 
        /home/sravisha/fuser/cNR/scripts/bedClip ${fe_path}.genome ${fe_path} ${fe_path}.clip
        LC_COLLATE=C sort -k1,1 -k2,2n ${fe_path}.clip > ${fe_path}.sort.clip
        /home/sravisha/fuser/cNR/scripts/bedGraphToBigWig ${fe_path}.sort.clip ${cl_path} ${sample}_FE.bw
        rm -f ${fe_path}.clip ${fe_path}.sort.clip

        bedtools slop -i ${llr_path} -g ${cl_path} -b 0 > ${llr_path}.genome 
        /home/sravisha/fuser/cNR/scripts/bedClip ${llr_path}.genome ${llr_path} ${llr_path}.clip
        LC_COLLATE=C sort -k1,1 -k2,2n ${llr_path}.clip > ${llr_path}.sort.clip
        /home/sravisha/fuser/cNR/scripts/bedGraphToBigWig ${llr_path}.sort.clip ${cl_path} ${sample}_logLR.bw
        rm -f ${llr_path}.clip ${llr_path}.sort.clip
        """
}

cutsite_path.into{ filter_path; plot_path }

process filterBigWig {
    echo false
    publishDir "$params.output.folder/finalBW", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "R/3.6.2-foss-2019b-fh1"

    input:
        set val(sample), path(fe_path), path(llr_path) from bwig_out
        each cut_site from filter_path
    
    output:
        tuple val(sample), path("${fe_path.getSimpleName()}_cutSite.bw"), path("${llr_path.getSimpleName()}_cutSite.bw") into fwig_out

    script:
        """
        #!/usr/bin/env Rscript
        library(rtracklayer)
        bigWigs <- import(\"${fe_path}\", format="BigWig")
        bedFile <- import(\"${cut_site}\", format="BED")
        seqlengths(bedFile) <- seqlengths(bigWigs)[names(seqlengths(bedFile))]
        overlap <- subsetByOverlaps(bigWigs, bedFile, type="any") 
        export(overlap, \"${fe_path.getSimpleName()}_cutSite.bw\", format="BigWig")
        bigWigs <- import(\"${llr_path}\", format="BigWig")
        bedFile <- import(\"${cut_site}\", format="BED")
        seqlengths(bedFile) <- seqlengths(bigWigs)[names(seqlengths(bedFile))]
        overlap <- subsetByOverlaps(bigWigs, bedFile, type="any") 
        export(overlap, \"${llr_path.getSimpleName()}_cutSite.bw\", format="BigWig")
        """
}

/*process plotPeaks {
    echo false
    publishDir "$params.output.folder/ngsPlots/${sample}", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    module "R/3.6.2-foss-2019b-fh1"

    input:
        set val(sample), path(bam_path), path(bai_path) from plot_out
        each cut_site from plot_path
        each region from plot_group
    
    output:
        tuple val(sample), path("${sample}*") into plot_output

    script:
        if( region == 'bed')
            """
            /home/sravisha/software/ngsplot-2.63/bin/ngs.plot.r -G hg38 -R bed -E ${cut_site} -C ${bam_path} -O ${sample}_bed 
            """
        else
            """
            /home/sravisha/software/ngsplot-2.63/bin/ngs.plot.r -G hg38 -R ${region} -C ${bam_path} -O ${sample}_${region}
            """

}*/

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