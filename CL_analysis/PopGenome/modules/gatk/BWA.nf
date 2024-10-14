process BWA {

    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }

    publishDir 'raw_vcfs', mode: 'copy', overwrite: true, pattern: '*g.vcf.gz'

    input:
    val(sample)

    output:
    tuple val(sample), file("${sample}.sam")

    script:
    """
    #!/bin/bash

    ${projectDir}/software/bwa/bwa mem -M -t 16 \
	${params.fasta} \
	${params.reads}${sample}/${sample}_220809_L004_R1.fastq.gz ${params.reads}${sample}/${sample}_220809_L004_R2.fastq.gz \
	> ${sample}.sam
	

    """
}

