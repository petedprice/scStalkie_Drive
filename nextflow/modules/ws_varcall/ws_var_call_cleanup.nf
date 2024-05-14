process ws_var_call_cleanup {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    input:
    tuple val(species), val(sample),  file("${sample}_crdata")

    output:
    tuple val(species), val(sample), file("${sample}_duplicates.bam")


    script:
    """
    #!/bin/bash
    
    gatk MarkDuplicates -h

    gatk MarkDuplicates \
	-I ${sample}_crdata/outs/possorted_genome_bam.bam \
	-O ${sample}_duplicates.bam \
	-M ${sample}_dup_metrics.txt \
	--REMOVE_DUPLICATES

    """
}

