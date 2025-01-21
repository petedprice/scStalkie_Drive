process seurat_markers_noX {

    label 'seurat'

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 340.GB }
    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: 'marker_seurat_noX.RData'

    input: 
    file("integrated_seurat_noX.RData")

    output:
    file("marker_seurat_noX.RData")

    script:
    """
    #!/bin/bash
    echo dog
    Rscript ${projectDir}/Rscripts/seurat/cell_type_ID.R \
	integrated_seurat_noX.RData \
	. \
	${params.celltype_markers}

    mv outdata/marker_seurat.RData marker_seurat_noX.RData	


    """
}


