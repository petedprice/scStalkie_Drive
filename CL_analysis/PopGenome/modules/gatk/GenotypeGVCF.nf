process GenotypeGVCF {

    label 'gatk'
    errorStrategy 'retry'

    tag {'genotype_' + contig + ctg_st + ctg_nd }

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    input:
    tuple file("genomics_db"), val(contig), val(ctg_st), val(ctg_nd)

    output:
    tuple file("${contig}_${ctg_st}_${ctg_nd}.vcf.gz"), val(contig)

    script:
    """
    #!/bin/bash
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
	-R ${params.fasta} \
	-V gendb://genomics_db \
	-O ${contig}_${ctg_st}_${ctg_nd}.vcf.gz

    """
}

