process snpeff_eff {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'snpeff'

    tag {'snpeff_eff_' + contig }

    //publishDir 'vcftools_filt_vcfs', mode: 'copy', overwrite: true, pattern: '*filt_fin.vcf.gz'

    input:
    tuple file("${contig}_filt.vcf.gz"), val(contig)

    output:
    tuple file("${contig}_synnon.vcf"), val(contig)

    
    script:
    """
    #!/bin/bash

    snpeff -eff ${params.snpeff_db} ${contig}_filt.vcf.gz > tmp.vcf
    grep -E '^#|synonymous_variant|missense_variant' tmp.vcf > ${contig}_synnon.vcf


    """
}
