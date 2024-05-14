process R_var_filtering {

    label 'tidyverse'

    cpus { 4 * task.attempt } 
    errorStrategy 'retry'
    maxRetries 6
    memory { 512.GB * task.attempt }


    publishDir 'filtered_var', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(species), val(sample),  file("${species}_${sample}_fin_snp_read.txt.gz")

    output:
    tuple val(sample), file("${sample}_snp_summarise.txt.gz")
    script:
    """
    #!/bin/bash

    Rscript ${projectDir}/Rscripts/var_filtering.R \
        ${species}_${sample}_fin_snp_read.txt.gz \
	$sample


    """
}
