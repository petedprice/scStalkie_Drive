process ws_var_call_HC {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }

    publishDir 'raw_vcfs', mode: 'copy', overwrite: true, pattern: '*g.vcf.gz'

    input:
    tuple val(species), val(sample), file(bam), val(ref)

    output:
    tuple val(species), val(sample), file(bam), file("${sample}.g.vcf.gz"), val(ref)

    script:
    """
    #!/bin/bash
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.fasta_dir}/${ref}.fna -I $bam -O ${sample}.g.vcf.gz -ERC GVCF
    """
}

