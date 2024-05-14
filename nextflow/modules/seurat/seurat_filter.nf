process seurat_filter {

    label 'seurat'

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 96.GB * task.attempt }

    publishDir 'initial_QC', mode: 'copy', overwrite: true, pattern: '*QC_plots'

    input: 
    tuple val(species), val(sample), file("${sample}_crdata"), val(ref), val(mt), val(nfs), val(mtr), val(gu)

    output:
    tuple val(sample), file("filtered_seurat.RData"), file("${sample}_initial_QC_plots"), val(nfs), val(mtr), val(gu)   
    
    script:
    """
    #!/bin/bash

    cat ${params.gtf_dir}/${ref}.gtf | grep $mt | grep gene_id > ss.gtf
    cut -f9 ss.gtf | cut -d ";" -f1 | cut -d " " -f2 | uniq | sed 's/"//g' > mt_genes.txt


    cat ${params.metadata} | grep ${sample} > metadata_ss.csv   

    Rscript ${projectDir}/Rscripts/seurat/seurat_filter.R \
	 ${projectDir} \
	 ${sample}_crdata/outs/filtered_feature_bc_matrix \
	 . \
	metadata_ss.csv \
	mt_genes.txt \
	$nfs \
	$mtr \
	$gu
    

    mv plots ${sample}_initial_QC_plots

    mv outdata/filtered_seurat.RData .


    """
}


