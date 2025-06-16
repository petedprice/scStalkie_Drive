process prank_phy {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'prank_convert_to_phy' + '_' + og }

    publishDir 'unmasked_allignments', mode: 'copy', overwrite: true, pattern: '*.phy'

    label 'prank'

    input:
    tuple file("${og}_codon.nogaps.fas"), val(og)

    output:
    tuple file("${og}_codon.nogaps.phy"), val(og)
 
    script:
    """
    #!/bin/bash
    ${baseDir}/software/prank/bin/prank -convert -d=${og}_codon.nogaps.fas -f=phylips -o=${og}_codon.nogaps
    """
}
