process seurat_markers {

    label 'seurat'

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 340.GB }
    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: 'marker_seurat.RData'

    input: 
    file("integrated_seurat.RData")

    output:
    file("marker_seurat.RData")

    script:
    """
    #!/bin/bash
    Rscript ${projectDir}/Rscripts/seurat/cell_type_ID.R \
	integrated_seurat.RData \
	. \
	${params.celltype_markers}

    mv outdata/marker_seurat.RData .	


    """
}


