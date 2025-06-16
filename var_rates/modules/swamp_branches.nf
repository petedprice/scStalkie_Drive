process mod0_paml {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'swamp_branches' + '_' + og }

    input:
    tuple file("${og}_swamp_analysis"), val(og)

    output:
    tuple file("${og}_swamp_analysis"), val(og), file('SWAMP_BRANCHES.txt')

    script:
    """
    #!/bin/bash

    Rscript ${baseDir}/scripts/branch_file.R ${og}_swamp_analysis/${og}_mod0.txt
    
    """
}
