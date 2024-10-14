process NCgene_BED {
    cpus = 4
    memory = '4 GB'
    time = '4h'

    label 'bedops'

    tag {'select_introns__' + contig }

    input:
    tuple val(vcf), val(contig), val(ctg_len)

    output:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_NC.bed")    
    
    script:
    """
    #!/bin/bash
    cat ${params.gff} | grep $contig | 
	awk -F'\t' '\$3 == "intron"' | gff2bed > ${contig}_NC.bed
    """
}
