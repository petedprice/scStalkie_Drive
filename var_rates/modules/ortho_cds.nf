process ortho_cds {
    cpus = 4
    memory = '64 GB'
    time = '8h'


    input:
    tuple file(proteins), file('Orthofinder_Results')

    output:
    tuple file('*_folder'), file('species_order.txt')


    script:
    """
    #!/bin/bash
    cat *cds_longest.fna > comp_longest.fna
    ls Orthofinder_Results/Single_Copy_Orthologue_Sequences | sed "s/.fa//" > sc_ogs.txt
    head -1 Orthofinder_Results/Phylogenetic_Hierarchical_Orthogroups/N0.tsv | cut -f4- > species_order.txt
    for c in \$(cat sc_ogs.txt )
    do
    mkdir \${c}_folder
    grep -w \$c Orthofinder_Results/Phylogenetic_Hierarchical_Orthogroups/N0.tsv | cut -f4- | tr '\\t' '\\n' | dos2unix > tmp_orthos_name.txt
    grep --no-group-separator -A 1 -F -f tmp_orthos_name.txt comp_longest.fna > \${c}_cds.fa
    mv \${c}_cds.fa \${c}_folder
    cp Orthofinder_Results/Single_Copy_Orthologue_Sequences/\${c}.fa \${c}_folder/\${c}_pro.fa    
    done
    """

}
