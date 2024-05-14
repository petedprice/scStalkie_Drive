process cluster_var_call_VF {

    label 'vcftools'
    errorStrategy 'retry'

    cpus { 2 * task.attempt }
    memory { 8.GB * task.attempt }

    publishDir 'filtered_vcfs_clusters', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:
    tuple val(sample), val(species), val(clus), file("${sample}_dup_NCR_clus${clus}.bam"), file("${sample}_clus${clus}.vcf.gz"), val(ref)

    output:
    tuple val(species), val(sample), val(clus), file("${sample}_dup_NCR_clus${clus}.bam"), file("${sample}_clus${clus}_filt.recode.vcf"), val(ref)
    script:
    """
    #!/bin/bash

    #gatk IndexFeatureFile -F ${sample}.g.vcf.gz

    vcftools \
	--gzvcf ${sample}_clus${clus}.vcf.gz \
	--out ${sample}_clus${clus}_filt \
	--minGQ 20 \
	--minDP 2 \
	--minQ 30 \
	--recode \
	--recode-INFO-all

    """
}

