process CollectWgsMetrics {

    label 'gatk'
    errorStrategy 'retry'

    tag {'haplotypecaller_' + sample + contig + ctg_st + '_' + ctg_nd }

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }
    time = '12h'  

    //publishDir 'raw_vcfs', mode: 'copy', overwrite: true, pattern: '*g.vcf.gz'

    input:
    tuple val(sample), file("${sample}_dup_NCR_RG.bam"), file("${sample}_dup_NCR_RG.bai"),  val(contig), val(ctg_st), val(ctg_nd)


    output:

    script:
    """
    #!/bin/bash

    gatk CollectWgsMetrics \
	-I ${sample}_dup_NCR_RG.bam \
	-O ${sample}_wgs_metrics.txt \
	-R ${params.fasta}

    """
}

