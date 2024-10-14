process MergeVcfs {

    label 'gatk'
    errorStrategy 'retry'

    tag {'mergevcfs_' + ctg  }

    cpus { 2 * task.attempt }
    errorStrategy 'retry'
    memory { 16.GB * task.attempt }

    input:
    tuple file(vcf), val(contig)

    output:
    tuple file("${contig}.vcf.gz"), val(contig)

    script:
    """
    #!/bin/bash

    ls -1 *vcf.gz > input_variant_files.list

    gatk MergeVcfs \
          -I input_variant_files.list \
          -O ${contig}.vcf.gz

    """
}

