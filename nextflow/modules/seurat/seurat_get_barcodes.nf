process seurat_get_barcodes {

    label 'seurat'

    cpus { 2 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 8.GB }

    input: 
    file("integrated_seurat.RData")

    output:
    file("*cluster.txt")

    script:
    """
    #!/bin/bash

    Rscript ${projectDir}/Rscripts/seurat/seurat_print_barcodes.R \
	integrated_seurat.RData

    """
}


