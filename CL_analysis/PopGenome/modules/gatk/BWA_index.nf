process BWA_index {

    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }

    publishDir 'raw_vcfs', mode: 'copy', overwrite: true, pattern: '*g.vcf.gz'

    input:

    output:

    script:
    """
    #!/bin/bash

    ${projectDir}/software/bwa/bwa index ${params.fasta}

    """
}

