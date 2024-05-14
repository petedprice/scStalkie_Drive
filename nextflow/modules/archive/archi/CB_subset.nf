process CB_subset {

    input:
    tuple val(species), val(sample), val(contig), file("subset.bam"), file("subset.bam.bai"), val(ref)

    output:
    tuple val(species), val(sample), val(contig), file("CB_subset.bam"), file("CB_subset.bam.bai"), env(ref_genome)    

    """
    #!/bin/bash
    
    ref_genome=${params.fasta_dir}/${ref}.fna
    
    #Filter bam for only barcoded reads
    samtools view -b -d CB -h -F 1024 -q 255 subset.bam > CB_subset.bam
    samtools index CB_subset.bam

    """
}
