process COUNT_SNPS {
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    maxRetries 8
    cpus = { 1 }
    memory = { 128.GB * task.attempt }
    time = {4.hour } //* task.attempt }

    tag {'PopGenome_' + contig + '_' + bedtype}

    publishDir 'PopGenome_Stats', mode: 'copy', overwrite: true, pattern: '*_PopGenome_Stats.csv'

    input:
    tuple val(contig), val(bedtype), file("subsettedAll_PH_${bedtype}_${contig}_PopGenom_Stats.csv"), file(vcf)

    output:
    tuple val(contig), val(bedtype), file("${contig}_${bedtpye}_No_SNPS.txdt") 

    script:
    """
    #!/bin/bash
    zcat $vcf | grep -v '#' | wc -l > ${contig}_${bedtype}_No_SNPS.txt    
    """
}
