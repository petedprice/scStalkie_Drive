process ws_var_call_splitNCR {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    publishDir 'dup_NCR_clean_bam', mode: 'copy', overwrite: true, pattern: '*dup_NCR.bam'

    input:
    tuple val(species), val(sample), file("${sample}_duplicates.bam"), val(ref)

    output:
    tuple val(species), val(sample), file("${sample}_dup_NCR.bam"), val(ref)

    script:
    """
    #!/bin/bash
    
    gatk SplitNCigarReads \
      -R ${params.fasta_dir}/${ref}.fna \
      -I ${sample}_duplicates.bam \
      -O ${sample}_dup_NCR.bam

    """
}

