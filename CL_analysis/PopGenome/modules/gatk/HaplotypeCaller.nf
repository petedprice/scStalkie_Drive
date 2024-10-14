process HaplotypeCaller {

    label 'gatk'
    errorStrategy 'retry'

    tag {'haplotypecaller_' + sample + contig + ctg_st + '_' + ctg_nd }

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }
    time = '4h'  

    //publishDir 'raw_vcfs', mode: 'copy', overwrite: true, pattern: '*g.vcf.gz'

    input:
    tuple val(sample), file("${sample}_dup_NCR_RG.bam"), file("${sample}_dup_NCR_RG.bai"),  val(contig), val(ctg_st), val(ctg_nd)

    output:
    tuple file("${sample}_${contig}_${ctg_st}_${ctg_nd}.g.vcf.gz"), val(contig), val(ctg_st), val(ctg_nd)

    script:
    """
    #!/bin/bash

    [[ $contig == "Chr_X" ]] && ploidy=1 || ploidy=2

    gatk --java-options "-Xmx110g" HaplotypeCaller \
	--native-pair-hmm-threads 6 \
	-R ${params.fasta} \
	-I ${sample}_dup_NCR_RG.bam \
	-O ${sample}_${contig}_${ctg_st}_${ctg_nd}.g.vcf.gz \
	-L ${contig}:${ctg_st}-${ctg_nd} \
	-ploidy \$ploidy \
	-ERC GVCF

    """
}

