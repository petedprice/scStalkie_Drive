process vcftools_filt {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'vcftools'

    tag {'vcftools_filt_' + contig }

    //publishDir 'vcftools_filt_vcfs', mode: 'copy', overwrite: true, pattern: '*filt_fin.vcf.gz'

    input:
    tuple file("${contig}_filt.vcf.gz"), val(contig)

    output:
    tuple file("${contig}_filt_fin.vcf.gz"), val(contig)

    
    script:
    """
    #!/bin/bash

    vcftools --gzvcf ${contig}_filt.vcf.gz --recode \
	--minDP 5 --minGQ 20 --max-missing 0.4 --max-alleles 2 --remove-indels \
	--stdout | gzip -c > ${contig}_filt_fin.vcf.gz

    """
}
