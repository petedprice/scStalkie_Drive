process alevin_fry_quant {
    cache 'lenient' 
    cpus = 16
    memory = '64 GB'
    time = '4h'


    tag {'Cellranger_count_' + '_' + species + '_' + sample }

    publishDir 'cellranger_out', mode: 'copy', overwrite: true, pattern: '*crdata'

    input:
    tuple val(species), file("${species}_index"), val(sample)

    output:
    //tuple val(species), val(sample), file("${sample}_crdata")

    script:
    """
    #!/bin/bash
    salmon alevin \
	-lISR --chromium \
	-1 <read1_files> \
	-2 <read2_files> \
	-o <alevin_odir> \
	-i ${species}_index \
	-p $tastk.cpus \
	--sketch


    """

}
