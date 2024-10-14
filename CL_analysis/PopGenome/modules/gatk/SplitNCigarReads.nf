process SplitNCigarReads {

    label 'gatk'
    errorStrategy 'retry'

    tag {'splitNCR_' + sample }

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    publishDir 'dup_NCR_clean_bam', mode: 'copy', overwrite: true, pattern: '*dup_NCR.bam'

    input:
    tuple val(sample), file("${sample}_marked_duplicates.bam")

    output:
    tuple val(sample), file("${sample}_dup_NCR.bam")

    script:
    """
    #!/bin/bash
    
    gatk SplitNCigarReads \
      -R ${params.fasta} \
      -I ${sample}_marked_duplicates.bam \
      -O ${sample}_dup_NCR.bam

    """
}

