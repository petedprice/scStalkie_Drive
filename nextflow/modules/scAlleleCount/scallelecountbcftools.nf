process ws_var_call_VF_stringent {

    label 'vcftools'
    errorStrategy 'retry'

    cpus { 4 * task.attempt }
    memory { 32.GB * task.attempt }

    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:
    tuple val(species), val(sample), file("${sample}_dup_NCR.bam"), file("${sample}.g.vcf.gz"), val(ref)

    output:
    tuple val(species), val(sample), file("${sample}_filt.recode.vcf"), val(ref), file("${sample}_dup_NCR.bam")

    script:
    """
    #!/bin/bash

    #gatk IndexFeatureFile -F ${sample}.g.vcf.gz

    vcftools \
	--gzvcf ${sample}.g.vcf.gz \
	--out ${sample}_filt \
	--minGQ 30 \
	--minDP 10 \
	--minQ 500 \
	--recode \
	--recode-INFO-all


    """
}

