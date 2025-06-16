process m1avsm2a {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'m1avsm2a' + '_' + og }

    publishDir 'pergene_m1am2a_summaries', mode: 'copy', overwrite: true, pattern: '*_paml_mod1a2a.txt'

    label 'R'

    input:
    tuple file(phy), val(og), file('tree_paml.txt')
    
    output:
    tuple file("*.txt"), env(sc)
 
    script:
    """
    #!/bin/bash
    cp ${baseDir}/data/PAML_CTLs/mod1a2a.ctl .

    sed -i 's/ALN/$phy/g' mod1a2a.ctl
    sed -i 's/TREE/tree_paml.txt/g' mod1a2a.ctl
    out="\$(basename /$phy .phy)"
    sc=\${out#"${og}_codon.nogaps_"}

    sed -i 's/OUT/paml_mod1a2a.txt/g' mod1a2a.ctl

    ${baseDir}/software/paml4.8/bin/codeml mod1a2a.ctl
    mv paml_mod1a2a.txt \${out}_paml_mod1a2a.txt


    """
}

