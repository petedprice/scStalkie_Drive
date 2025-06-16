process comp_paml_models {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'comp_paml_models' + '_' + sc }

    publishDir 'model_summaries', mode: 'copy', overwrite: true, pattern: '*_model_comparison.txt'

    label 'R'

    input:
    tuple file(mod1a2a), val(sc)
    
    output:
    file("*${sc}_*")
 
    script:
    """
    #!/bin/bash

    Rscript ${baseDir}/scripts/comp_paml_models.R . 2 ${sc}_model_comparison.txt

    """
}

