process FD_BED {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    tag {'select_sites_FourD_' + contig }

    input:
    tuple val(vcf), val(contig), val(ctg_len)

    output:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_FD.bed")

    script:
    """
    #!/bin/bash
    cat ${params.degen} | grep $contig | 
	awk -F'\t' '\$5 == "4"' > ${contig}_FD.bed
    """
}
