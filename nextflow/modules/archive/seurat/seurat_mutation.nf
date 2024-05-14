process seurat_mutation {

    label 'seurat'

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 360.GB }

    publishDir 'mut_metadata', mode: 'copy', overwrite: true, pattern: '*_mutant_metadata.csv'
    publishDir 'mut_metadata', mode: 'copy', overwrite: true, pattern: '*_germ_somatic.txt'
    publishDir 'seurat_RData', mode: 'copy', overwrite: true, pattern: '*_marker_seurat.RData'


    input: 
    tuple val(samples), val(sex), val(stage), val(species), file(mutation_data), file("marker_seurat.RData")

    output:
    tuple val(samples), val(sex), val(stage), val(species), file("*_mutant_metadata.csv"), file("*_marker_seurat.RData"), file("*_germ_somatic.txt")

    script:
    """
    #!/bin/bash

    cat ${params.metadata} | grep $species | grep ,$sex, | grep $stage > metadata_ss.csv
    
    echo blobs
  
    echo $samples | sed 's/[],[]//g' > samples.txt
    for samp in \$(cat samples.txt)
	do 
	echo \$samp > tmp.txt
	cat metadata_ss.csv | grep \$samp, > md_ss.csv
	cat md_ss.csv
	    Rscript ${projectDir}/Rscripts/seurat/4.seurat_mutation.R \
		marker_seurat.RData \
		.\
		tmp.txt \
		md_ss.csv
    
	    #mv outdata/seurat_mutant.RData \${s}_${species}_${sex}_${stage}_seurat_mutant.RData 
	done
    mv marker_seurat.RData ${species}_${sex}_${stage}_marker_seurat.RData

    """
}


