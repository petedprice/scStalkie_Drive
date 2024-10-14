process tabix {
    cpus = 8
    memory = '64 GB'
    time = '4h'

    label 'vcftools'

    tag {'tabix_' + contig}

    input:
    tuple file("${contig}_filt_fin.vcf.gz"), val(contig)

    output:
    tuple val(contig), file("${contig}_filt_tbx.vcf.gz"), file("${contig}_filt_tbx.vcf.gz.tbi")
    
    script:
    """
    #!/bin/bash

    cp ${contig}_filt_fin.vcf.gz ./${contig}_filt_tbx.vcf.gz
    gunzip ${contig}_filt_tbx.vcf.gz
    bgzip ${contig}_filt_tbx.vcf
    tabix -p vcf ${contig}_filt_tbx.vcf.gz

    """
}
