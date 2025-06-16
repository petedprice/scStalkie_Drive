process remove_Ns {
    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 1.GB * task.attempt }

    tag {'remove_Ns' + '_' + og }

    publishDir 'NoNs_allignments', mode: 'copy', overwrite: true, pattern: '*NoNs.phy'

    label 'R'

    input:
    tuple file(phy), val(og)

    output:
    tuple file('*NoNs.phy'), val(og)

    script:
    """
    #!/bin/bash
    Rscript ${baseDir}/scripts/remove_Ns.R $phy 150



    """
}
