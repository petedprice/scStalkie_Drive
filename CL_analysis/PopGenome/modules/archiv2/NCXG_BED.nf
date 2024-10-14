process NCXG_BED {
    cpus = 4
    memory = '4 GB'
    time = '4h'

    label 'bedops'

    tag {'select_introns__' + contig }

    input:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_NC.bed"), file("${contig}_Xgene.bed")

    output:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_NCXG.bed")    
    
    script:
    """
    #!/bin/bash
    bedops --difference ${contig}_Xgene.bed ${contig}_NC.bed > ${contig}_NCXG.bed

    """
}
