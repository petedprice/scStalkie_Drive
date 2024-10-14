process pre_HaplotypeCaller {

    label 'gatk'
    errorStrategy 'retry'

    tag {'haplotypecaller_' + sample + contig + ctg_st + '_' + ctg_nd }

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }
    time = '12h'  

    //publishDir 'raw_vcfs', mode: 'copy', overwrite: true, pattern: '*g.vcf.gz'

    input:
    tuple val(sample), file("${sample}_dup_NCR.bam"), val(contig), val(ctg_st), val(ctg_nd)

    output:
    tuple val(sample), file("${sample}_dup_NCR_RG.bam"), file("${sample}_dup_NCR_RG.bai"),  val(contig), val(ctg_st), val(ctg_nd)

    script:
    """
    #!/bin/bash

    gatk AddOrReplaceReadGroups \
	-I ${sample}_dup_NCR.bam \
	-O ${sample}_dup_NCR_RG.bam \
	-SM $sample \
	-LB TBC \
	-PL TBC \
	-PU TBC

    gatk BuildBamIndex \
	-I ${sample}_dup_NCR_RG.bam

    """
}

