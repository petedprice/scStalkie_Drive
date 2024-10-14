process MKT_vcftools {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'vcftools'

    tag {'vcftools_MKT_filt_' + contig }

    //publishDir 'vcftools_filt_vcfs', mode: 'copy', overwrite: true, pattern: '*filt_fin.vcf.gz'

    input:
    tuple file("${contig}_synnon.vcf"), val(contig)

    output:
    tuple file("${contig}_synnon_maf10.vcf"), val(contig)

    
    script:
    """
    #!/bin/bash

    vcftools \
	--vcf ${contig}_synnon.vcf \
	--maf 0.10 --recode-INFO-all \
	--recode --out ${contig}_synnon_maf10


    cat ${params.sample_info} | grep ST | cut -f1 -d ',' > ST_samples.txt
    cat	${params.sample_info} |	grep SR	| cut -f1 -d ',' > SR_samples.txt

    vcftools --vcf ${contig}_synnon_maf10.vcf --keep ST.txt --freq --out ${contig}_ST_freq
    vcftools --vcf ${contig}_synnon_maf10.vcf --keep SR.txt --freq --out ${contig}_SR_freq


    """
}
