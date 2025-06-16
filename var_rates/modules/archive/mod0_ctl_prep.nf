process mod0_ctl_prep {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    label 'paml'

    input:

    output:
    file('dog.txt')

    script:
    """
    #!/bin/bash
    touch dog.txt    

    """
}
