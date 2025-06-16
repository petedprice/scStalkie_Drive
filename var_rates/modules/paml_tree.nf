process paml_tree {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    label 'R'

    input:
    file(tree)

    output:


    file('tree_paml.txt')
 
    script:
    """

    #!/bin/bash

    cp ${params.tree} tree.txt

    Rscript ${baseDir}/scripts/paml_tree.R tree.txt







    """
}

