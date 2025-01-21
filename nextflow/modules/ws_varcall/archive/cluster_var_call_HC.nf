process cluster_var_call_HC {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 2 * task.attempt }
    errorStrategy 'retry'
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample), file(cluster), val(clus), val(species), file("${sample}_dup_NCR_clus${clus}.bam"), val(ref)

    output:
    tuple val(sample), val(species), val(clus), file("${sample}_dup_NCR_clus${clus}.bam"), file("${sample}_clus${clus}.vcf.gz"), val(ref)

    script:
    """
    #!/bin/bash
    samtools index ${sample}_dup_NCR_clus${clus}.bam
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.fasta_dir}/${ref}.fna -I ${sample}_dup_NCR_clus${clus}.bam -O ${sample}_clus${clus}.vcf.gz -ERC NONE
    """
}

