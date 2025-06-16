process cellranger_mkref {

    cpus = 16
    memory = '64 GB'
    time = '4h'

    input:
    tuple val(species), val(ref)
    output:
    tuple val(species), file("${species}_cellranger_reference")
    script:
    """
    #!/bin/bash

    $params.cellranger mkgtf \
	${params.gtf_dir}/${ref}.gtf \
	cr_${ref}.gtf

    $params.cellranger mkref \
	--genome=${species}_cellranger_reference \
	--fasta=${params.fasta_dir}/${ref}.fna \
        --genes=cr_${ref}.gtf


    """
}
