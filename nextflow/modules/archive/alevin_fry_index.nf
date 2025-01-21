process alevin_fry_index {

    cpus = 4
    memory = '16 GB'
    time = '4h'

    //conda "bioconda::alevin-fry bioconda::salmon" 

    input:
    tuple val(species), val(ref)
    output:
    tuple val(species), file("${species}_index")
    script:
    """

    #!/bin/bash
    module load miniconda
    source activate alevin
    salmon index -t ${params.alevin_genome}/${ref}_cds.fa -i ${species}_index


    """
}
