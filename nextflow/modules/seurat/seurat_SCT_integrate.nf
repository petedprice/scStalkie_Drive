process seurat_SCT_integrate {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 84.GB * task.attempt }    
    
    label 'seurat'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple file("doublet_seurat.RData"), val(nfs), val(mtr), val(gu)

    output:
    tuple file("integrated_seurat_nf${nfs}_mtr${mtr}_gu${gu}.RData"), val(nfs), val(mtr), val(gu)

    script:
    """
    #!/bin/bash
    echo ${task.memory}
    Rscript ${projectDir}/Rscripts/seurat/SCT_integrate.R \
	.

 	
    mv outdata/integrated_seurat.RData integrated_seurat_nf${nfs}_mtr${mtr}_gu${gu}.RData

    """
}


