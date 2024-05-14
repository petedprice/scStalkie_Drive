process seurat_SCT_integrate {

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 340.GB }
    
    label 'seurat'

    //publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(samples), val(sex), val(stage), val(species), val(ref), file("doublet_seurat.RData")

    output:
    tuple val(samples), val(sex), val(stage), val(species), val(ref), file("integrated_seurat.RData")

    script:
    """
    #!/bin/bash
    echo ${task.memory}
    Rscript ${projectDir}/Rscripts/seurat/subsetted_funcs/SCT_integrate.R \
	.

 	
    mv outdata/integrated_seurat.RData .

    """
}


