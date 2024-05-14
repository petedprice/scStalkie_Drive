process mut_compiled {

    publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*mutation_data.txt.gz'

    input:
    tuple val(sample), file(snps)

    output:
    tuple val(sample), file("${sample}_mutation_data.txt.gz")
    
    script:
    """
    #!/bin/bash
    cat *.txt.gz > ${sample}_mutation_data.txt.gz
    """
}
