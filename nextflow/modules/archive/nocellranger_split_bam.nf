process nocellranger_split_bam {
 

    label 'samtoolsetc' 
    

    input:
    tuple val(species), val(sample), file("contig.txt")

    output:
    tuple val(species), val(sample), env(contig), file("subset.bam"), file("subset.bam.bai")

    script: 
    """
    #module load apps/SAMtools/1.7/gcc-4.9.4
    
    contig=\$(cat contig.txt) 
    echo \$contig
    samtools view -bh ${params.cellranger_data}/${sample}/outs/possorted_genome_bam.bam \$contig > subset.bam
    samtools index subset.bam    
    """
}
