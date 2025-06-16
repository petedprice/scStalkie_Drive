process format_fastas {
    cpus = 1
    memory = '1 GB'
    time = '1h'

    label 'R'

    input:
    tuple file('seq_folder'), file('species_order.txt')

    output:
    tuple file('*_fastas'), env(og)

    script:
    """
    #!/bin/bash
    cp seq_folder/*cds.fa cds.fa
    cp seq_folder/*pro.fa pro.fa

    Rscript ${baseDir}/scripts/rename_fasta.R \
	species_order.txt \
	cds.fa \
	cds_species.fa
	

    Rscript ${baseDir}/scripts/rename_fasta.R \
        species_order.txt \
        pro.fa \
        pro_species.fa

    
    ogt1==\$(ls seq_folder/*pro.fa)
    ogt2=\$(echo \${ogt1##*/})
    og=\${ogt2%_pro.fa}    

    mkdir \${og}_fastas

    mv *.fa \${og}_fastas
    mv species_order.txt \${og}_fastas

    """
}
