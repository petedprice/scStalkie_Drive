process sc_var_call {

    tag {'varcall_' + sample}

    label 'samtoolsetc'

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 7
    memory { 104.GB * task.attempt }


    //publishDir 'mut_id', mode: 'copy', overwrite: true, pattern: '*fin_snp*'


    input:
    tuple val(species), val(sample), file("${sample}_dup_NCR.bam"), file("${sample}_filt.recode.vcf"), val(ref)

    output:
    tuple val(species), val(sample), file("${species}_${sample}_fin_snp_read.txt.gz")

    """
    #!/bin/bash
    cp ${sample}_filt.recode.vcf snps.vcf
    gzip snps.vcf
    ref_genome=${params.fasta_dir}/${ref}.fna

    #Filter bam for only barcoded reads
    samtools index ${sample}_dup_NCR.bam
    samtools view -b -d CB -h -q 30 ${sample}_dup_NCR.bam > CB_subset.bam
    samtools index CB_subset.bam

    #RUN SAM2TSV
    java -jar /opt/jvarkit/dist/jvarkit.jar sam2tsv -R \$ref_genome --regions snps.vcf.gz -N -o sam2tsv.tsv.gz CB_subset.bam

    #GREP ALL SITES from bcftools in samttsv output to subset it
    zcat snps.vcf.gz | cut -f2 | egrep -v "^#"  | egrep -v POS > snps.txt
    awk 'NR==FNR{snps[\$1]; next} \$8 in snps' snps.txt <(zcat sam2tsv.tsv.gz) > snps_sam2tsv.tsv
    cat snps_sam2tsv.tsv | cut -f1 | uniq > uniq_reads.txt

    #You now have read names supporting each variant in the file
    #You now need to associate the read names with the barcode.

    #Subset bam to just those reads
    samtools view CB_subset.bam | grep -Ff uniq_reads.txt > urCB_subset.sam

    # Extract barcodes from read names and pair with read name
    grep -o -E '.CB:Z.{0,19}' urCB_subset.sam | cut -f2 | cut -c 6- > barcodes.txt
    cut -f1 urCB_subset.sam > reads.txt
    paste reads.txt barcodes.txt  > reads_barcodes.txt

    #Merge reads_barcodes.txt with your snps_sam2tsv.tsv.gz file
    join reads_barcodes.txt snps_sam2tsv.tsv | gzip > ${species}_${sample}_fin_snp_read.txt.gz

    """
}
