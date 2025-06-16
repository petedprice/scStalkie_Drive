process swamp {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    label 'python'

    input:
    tuple file("${og}_swamp_analysis"), val(og), val(

    output:

    script:
    """
    #!/bin/bash

    python2.7 ${baseDir}/software/SWAMP/SWAMP.py \
	-i ${og}_swamp_analysis \
	-b ${baseDir}/data/SWAMP_BRANCHES.txt \
	-t 7 \
	-w 15 \
	-m 100 \
	> ${og}_swamp_t7w15.txt

	


    """
}
