process GenomicsDBImport {

    label 'gatk'
    errorStrategy 'retry'
    tag {'genomicsdb_' + contig + ctg_st + ctg_nd }


    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    input:
    tuple file(vcf), val(contig), val(ctg_st), val(ctg_nd)

    output:
    tuple file("genomics_db"), val(contig), val(ctg_st), val(ctg_nd)

    script:
    """
    #!/bin/bash

    for file in *.g.vcf.gz; do
	gatk IndexFeatureFile -F \$file
	sample=\$(echo "\$file" | cut -d'_' -f1)
	contig=\$(echo "\$file" | cut -d'_' -f2) 
	echo -e "\${sample}\t\${file}" >> sample-name-map
    done


    gatk GenomicsDBImport \
	--sample-name-map sample-name-map \
	--genomicsdb-workspace-path genomics_db \
	--tmp-dir . \
	-L $contig    

    """
}

