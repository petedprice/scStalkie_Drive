process degenotate {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    input:

    output:
    
    script:
    """
    #!/bin/bash
    conda activate degenotate
    degenotate.py \
	  -a ${prams.gff} \
	  -g ${params.fasta} \
	  -o degenotate_data \
	  --overwrite
    conda deactivate

    """
}
