process bcftools_mut {

    input:
    tuple val(species), val(sample), val(contig), file("CB_subset.bam"), file("CB_subset.bam.bai"), val(ref_genome)

    output:
    tuple val(species), val(sample), val(contig), file("CB_subset.bam"), file("CB_subset.bam.bai"), val(ref_genome), file("snps.vcf.gz")

    """
    #!/bin/bash
    
    #Run BCFTOOLS
    bcftools mpileup --threads 1 -f $ref_genome CB_subset.bam | bcftools call --threads 1 -mv -Ob | bcftools view --types snps -i 'GT="het" &  INFO/DP >= 4 & (DP4[0]+DP4[1])>1 & (DP4[2]+DP4[3])>1' | gzip > snps.vcf.gz

    """
}
