process split_bam_cluster {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 2 * task.attempt }
    errorStrategy 'retry'
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample), file(cluster), val(clus), val(species), file("${sample}_dup_NCR.bam"), val(ref)

    output:
    tuple val(sample), file(cluster), val(clus), val(species), file("${sample}_dup_NCR_clus${clus}.bam"), val(ref)


    script:
    """
    #!/bin/bash

    samtools index ${sample}_dup_NCR.bam

    ${projectDir}/software/subset-bam_linux \
	--bam ${sample}_dup_NCR.bam \
	--cell-barcodes $cluster \
	--out-bam ${sample}_dup_NCR_clus${clus}.bam
    
    #samtools index ${sample}_dup_NCR_clus${clus}.bam

    #gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.fasta_dir}/${ref}.fna -I ${sample}_dup_NCR_clus${clus}.bam -O ${sample}.vcf.gz -ERC NONE
     
   
    """
}

