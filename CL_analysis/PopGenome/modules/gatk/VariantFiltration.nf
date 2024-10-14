process VariantFiltration {

    label 'gatk'
    errorStrategy 'retry'

    tag {'variantfilter_gatk_' + ctg }

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }
    time = '12h'  

    //publishDir 'gatk_filt_vcfs', mode: 'copy', overwrite: true, pattern: '*filt.vcf.gz'

    input:
    tuple file("${contig}.vcf.gz"), val(contig)

    output:
    tuple file("${contig}_filt.vcf.gz"), val(contig)

    script:
    """
    #!/bin/bash

    gatk IndexFeatureFile -F ${contig}.vcf.gz

    gatk SelectVariants \
	-V ${contig}.vcf.gz \
	-select-type SNP \
	-O snps.vcf.gz


    gatk VariantFiltration \
	-R ${params.fasta} \
	-V snps.vcf.gz \
	-O snps2.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"


    gatk SelectVariants \
	-V snps2.vcf.gz\
	--exclude-filtered \
	-O ${contig}_filt.vcf.gz

    """
}

