process sam2tsv {

    input:
    tuple val(species), val(sample), val(contig), file("CB_subset.bam"), file("CB_subset.bam.bai"), val(ref_genome), file("snps.vcf.gz")

    output:
    tuple val(species), val(sample), val(contig), file("CB_subset.bam"), file("CB_subset.bam.bai"), val(ref_genome)

    """
    #!/bin/bash
     
    #RUN SAM2TSV
    java -jar jvarkit.jar sam2tsv -R $ref_genome --regions snps.vcf.gz -N -o sam2tsv.tsv.gz CB_subset.bam

    #GREP ALL SITES from bcftools in samttsv output to subset it
    zcat snps.vcf.gz | cut -f2 | egrep -v "^#"  | egrep -v POS > snps.txt
    zgrep -Ff snps.txt sam2tsv.tsv.gz  >  snps_sam2tsv.tsv
    cat snps_sam2tsv.tsv | cut -f1 | uniq > uniq_reads.txt


    """
}
