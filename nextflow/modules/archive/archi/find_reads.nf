process mut_id {

    input:
    tuple val(species), val(sample), val(contig), file("subset.bam"), file("subset.bam.bai"), val(ref)

    output:
    tuple val(species), val(sample), val(contig), file("${species}_${sample}_${contig}_fin_snp_read.txt.gz")
    

    """
    #!/bin/bash
    #You now have read names supporting each variant in the file
    #You now need to associate the read names with the barcode.

    #Subset bam to just those reads
    samtools view CB_subset.bam | grep -Ff uniq_reads.txt > urCB_subset.sam

    # Extract barcodes from read names and pair with read name
    grep -o -E '.CB:Z.{0,19}' urCB_subset.sam | cut -f2 | cut -c 6- > barcodes.txt
    cut -f1 urCB_subset.sam > reads.txt
    paste reads.txt barcodes.txt  > reads_barcodes.txt

    #Merge reads_barcodes.txt with your snps_sam2tsv.tsv.gz file
    join reads_barcodes.txt snps_sam2tsv.tsv | gzip > ${species}_${sample}_${contig}_fin_snp_read.txt.gz

    """
}
