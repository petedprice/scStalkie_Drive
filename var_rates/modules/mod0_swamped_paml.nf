process mod0_swamped_paml {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'mod0' + '_' + og }

    publishDir 'check_trees', mode: 'copy', overwrite: true, pattern: '*SWAMP_BRANCHES.txt'

    label 'R'

    input:
    tuple file(phy), val(og)

    output:

    script:
    """
    #!/bin/bash

    cp ${baseDir}/data/PAML_CTLs/mod0.ctl .

    sed -i 's/ALN/${phy}/g' mod0.ctl
    sed -i 's/TREE/${params.tree}/g' mod0.ctl
    sed -i 's/OUT/${og}_mod0.txt/g' mod0.ctl

    ${baseDir}/software/paml4.8/bin/codeml mod0.ctl


    """
}
