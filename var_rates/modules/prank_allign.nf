process prank_allign {

    cpus { 2 * task.attempt }
    errorStrategy { task.exitStatus in 100..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'prank_allign' + '_' + og }

    label 'prank'

    input:
    tuple file("${og}_fastas"), val(og), file(tree)

    output:
    tuple file("${og}_codon.best.fas"), val(og)

    script:
    """
    #!/bin/bash

    cp ${og}_fastas/* .
    ulimit -c unlimited
    ${baseDir}/software/prank/bin/prank -d=cds_species.fa -t=$tree -codon -F -o=${og}_codon -once



    """
}
