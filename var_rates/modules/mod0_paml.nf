process mod0_paml {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'mod0' + '_' + og }

    publishDir 'check_trees', mode: 'copy', overwrite: true, pattern: '*SWAMP_BRANCHES.txt'

    label 'R'

    input:
    tuple file("${og}_codon.nogaps.phy"), val(og), file('tree_paml.txt')

    output:
    tuple file("${og}_swamp_analysis"), val(og), file("${og}_SWAMP_BRANCHES.txt")

    script:
    """
    #!/bin/bash

    cp ${baseDir}/data/PAML_CTLs/mod0.ctl .

    sed -i 's/ALN/${og}_codon.nogaps.phy/g' mod0.ctl
    sed -i 's/TREE/tree_paml.txt/g' mod0.ctl
    sed -i 's/OUT/${og}_mod0.txt/g' mod0.ctl

    ${baseDir}/software/paml4.8/bin/codeml mod0.ctl

    mkdir ${og}_swamp_analysis
    cp 2NG* ${og}_mod0.txt ${og}_codon.nogaps.phy rst rst1 rub lnf 4fold.nuc ${og}_swamp_analysis

    Rscript ${baseDir}/scripts/branch_file.R tree_paml.txt ${og}_mod0.txt 
    mv SWAMP_BRANCHES.txt ${og}_SWAMP_BRANCHES.txt

    """
}
