process scAlleleCount {

    errorStrategy 'retry'

    cpus { 4 * task.attempt }
    memory { 32.GB * task.attempt }

    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:
    tuple val(species), val(sample), file("${sample}_filt.recode.vcf"), val(ref), file("${sample}_dup_NCR.bam")

    output:

    script:
    """
    #!/bin/bash


    """
}

