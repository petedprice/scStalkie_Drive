process seurat_doublet {

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 340.GB }
    
    label 'seurat'

    //publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(sample), val(sex), val(stage), val(species), val(ref), file("filtered_seurat.RData")

    output:
    tuple val(sample), val(sex), val(stage), val(species), val(ref), file("doublet_seurat.RData")

    script:
    """
    #!/bin/bash
    echo ${task.memory}
    Rscript ${projectDir}/Rscripts/seurat/subsetted_funcs/doublet_finder.R \
	filtered_seurat.RData \
	. \
	${task.cpus} \
	${params.cellcycle_markers} \
	TRUE \

 	
    mv outdata/doublet_seurat.RData .

    """
}


