process branch_model {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'branch_models' + '_' + og }

    publishDir 'pergene_branch_model_summaries', mode: 'copy', overwrite: true, pattern: '?????'

    label 'R'

    input:
    tuple file(phy), val(og), file('paml_branch_trees')
    
    output:
    tuple file("*_branch_models"), env(sc)
 
    script:
    """
    #!/bin/bash
    cp paml_branch_trees/* .


    for tree in *paml_branch.txt
    do
    out=\${tree%_paml_branch.txt}_\$(basename /$phy .phy)

    cp ${baseDir}/data/PAML_CTLs/branch.ctl \$out.ctl
    cp ${baseDir}/data/PAML_CTLs/branch_null.ctl \${out}_null.ctl

    sed -i 's/ALN/$phy/g' \$out.ctl
    sed -i "s/TREE/\$(echo \$tree)/g" \${out}.ctl
    #sed -i 's/TREE/\$tree/g' \$out.ctl
    sed -i 's/OUT/paml_branch.txt/g' \$out.ctl

    ${baseDir}/software/paml4.8/bin/codeml \$out.ctl
    mv paml_branch.txt \${out}_paml_branch.txt

    sed -i 's/ALN/$phy/g' \${out}_null.ctl
    sed -i "s/TREE/\$(echo \$tree)/g" \${out}_null.ctl
    #sed -i 's/TREE/\$tree/g' \${out}_null.ctl
    sed -i 's/OUT/paml_branch_null.txt/g' \${out}_null.ctl

    ${baseDir}/software/paml4.8/bin/codeml \${out}_null.ctl
    mv paml_branch_null.txt \${out}_paml_branch_null.txt

    rm \$tree
    done

    out="\$(basename /$phy .phy)"
    sc=\${out#"${og}_codon.nogaps_"}

    mkdir \${out}_branch_models
    cp *paml_branch_null.txt \${out}_branch_models 
    cp *paml_branch.txt \${out}_branch_models


    """
}

