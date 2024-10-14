process tabix {
    cpus = 8
    memory = '64 GB'
    time = '4h'

    label 'vcftools'

    tag {'tabix_' + contig}

    input:
    tuple file("${contig}_filt_fin.vcf.gz"), val(contig)

    output:
    tuple val(contig), file("${contig}_filt_fin.vcf.gz"), file("${contig}_filt_fin.vcf.gz.tbi")
    
    script:
    """
    #!/bin/bash

    cp ${contig}_filt_fin.vcf.gz tmp.vcf.gz
    gunzip tmp.vcf.gz
    bgzip tmp.vcf
    tabix -p vcf tmp.vcf.gz
    mv tmp.vcf.gz.tbi ${contig}_filt_fin.vcf.gz.tbi


    """
}
