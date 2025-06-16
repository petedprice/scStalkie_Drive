process paml_branch_tree {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    label 'R'

    input:
    tuple file(tree), file(branch_tree)
    
    output:
    file('paml_branch_trees')
 
    script:
    """
    #!/bin/bash    
    cp $tree tree.txt
    cp -r $branch_tree branch_tree
    Rscript ${baseDir}/scripts/paml_branch_model_trees.R branch_tree
    """
}

