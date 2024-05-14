process cellranger_count {
    cache 'lenient' 
    cpus = 16
    memory = '64 GB'
    time = '4h'


    tag {'Cellranger_count_' + '_' + species + '_' + sample }

    //publishDir 'cellranger_out', mode: 'copy', overwrite: true, pattern: '*crdata'

    input:
    tuple val(species), file("${species}_cellranger_reference"), val(sample)

    output:
    tuple val(species), val(sample), file("${sample}_crdata")

    script:
    """
    #!/bin/bash
    $params.cellranger count \
	--id=$sample \
   	--fastqs=${params.read_dir}/${sample} \
   	--transcriptome=${species}_cellranger_reference \
	--localcores=${task.cpus} \
	--localmem=${task.memory.giga}
    
    mv $sample ${sample}_crdata

    """

}
