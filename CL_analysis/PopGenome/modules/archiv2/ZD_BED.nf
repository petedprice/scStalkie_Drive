process ZD_BED {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    tag {'select_sites_ZeroD_' + contig }

    input:
    tuple val(vcf), val(contig), val(ctg_len)

    output:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_ZD.bed")    
    
    script:
    """
    #!/bin/bash
    cat ${params.degen} | grep $contig | 
	awk -F'\t' '\$5 == "0"' > ${contig}_ZD.bed
    """
}
