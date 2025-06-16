process remove_gaps {
    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 1.GB * task.attempt }

    tag {'remove_gaps' + '_' + og }

    label 'R'

    input:
    tuple file("${og}_codon.best.fas"),	val(og)

    output:
    tuple file("${og}_codon.nogaps.fas"), val(og), optional: true

    script:
    """
    #!/bin/bash
    Rscript ${baseDir}/scripts/remove_gaps.R ${og}_codon.best.fas 150 ${og}_codon.nogaps.fas
    """
}
