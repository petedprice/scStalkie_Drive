process R_mut_filtering {

    label 'tidyverse'

    cpus { 5 * task.attempt } 
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 68.GB * task.attempt }


    //publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(species), val(sample), val(contig), file("${species}_${sample}_${contig}_fin_snp_read.txt.gz")

    output:
    tuple val(sample), file("${sample}_${contig}_snp_summarise.txt.gz"), optional: true    
    script:
    """
    #!/bin/bash
    echo running again agaain
    Rscript ${projectDir}/Rscripts/mut_filtering.R \
	${species}_${sample}_${contig}_fin_snp_read.txt.gz \
	$sample \
	$contig 
    echo printing

    """
}
