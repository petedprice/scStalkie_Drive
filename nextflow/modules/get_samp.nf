process get_samp {
    cpus = 1
    memory = '4 GB'
    time = '4h'


    input:
    file(cluster)
    output:
    tuple file(cluster), env(samp), env(clus)


    script:
    """
    #!/bin/bash
    cat $cluster
    tmp=\$(echo $cluster)
    delim="_"
    samp="\${tmp%%\$delim*}"
    clustmp="\${tmp#*\$delim}"
    clus="\${clustmp%cluster*}"
    """
}
