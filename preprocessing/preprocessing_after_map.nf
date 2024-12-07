#!/usr/bin/env nextflow
/*
* AUTHOR: Georgette Femerling 
* VERSION: 1.0
*/

/* How to run:
* nextflow run preprocessing_after_map.nf --ref /path/to/Reference/*.fasta --bams /path/to/merged_bamfiles/*.bam
* or with config to change default running parameters:
* nextflow run preprocessing_after_map.nf -c preprocessing_after_map.config
*/

nextflow.enable.dsl=2

params.ref = "/path/to/Reference/*.fasta" // absolute path to reference in fasta format
params.bams = "/path/to/merged_bamfiles/*.bam" // absolute path to Merged Bamfiles

process {
    executor = 'slurm'
    queue = 'hologenomics'
    clusterOptions = '--mem-per-cpu=40G'
    time = '10:00:00'
    cpus = 1
}

executor {
    name = 'slurm'
    queueSize = 13
}

process prepareReference {
    module 'htslib:samtools'

    input:
    path ref

    output:
    path "${ref}.fai", emit: fai
    path "${ref}.dict", emit: dict
    publishDir "${ref.parent}", mode: 'copy'

    script:
    """
    samtools faidx $ref
    samtools dict $ref > ${ref.baseName}.dict
    """
}

process markDuplicates {
    module 'java:picard'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.merged.dedup.bam"), path("${sample}.merged.dedup.bam.bai")
    publishDir "Bamfiles/Deduped", mode: 'copy'

    script:
    """
    picard MarkDuplicates INPUT=$bam \
        OUTPUT=${sample}.merged.dedup.bam \
        METRICS_FILE=${sample}.metrics \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true \
        REMOVE_SEQUENCING_DUPLICATES=true

    samtools index ${sample}.merged.dedup.bam
    """
}

process realignIndels {
    module 'GATK/v3.8.1'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref
    path ref_fai
    path ref_dict

    output:
    tuple val(sample), path("${sample}.merged.dedup.realigned.bam")
    publishDir "Bamfiles/Realigned", mode: 'copy'

    script:
    """
    gatk -T RealignerTargetCreator \
        -R $ref -I $bam \
        -o ${sample}_realignertargetcreator.intervals

    gatk -T IndelRealigner \
        -R $ref -I $bam \
        -targetIntervals ${sample}_realignertargetcreator.intervals \
        -o ${sample}.merged.dedup.realigned.bam
    """
}

workflow {
    reference = channel.fromPath(params.ref)
    bams = channel.fromPath(params.bams).map { file -> tuple(file.simpleName, file) }

    prepareReference(reference)
    markDuplicates(bams)
    realignIndels(markDuplicates.out, reference, prepareReference.out.fai, prepareReference.out.dict)
}
