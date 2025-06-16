process prank {
    cpus = 8
    memory = '32 GB'
    time = '4h'

    label 'prank'

    input:
    tuple file("seq_folder"), file('species_order.txt')


    output:


    script:
    """
    #!/bin/bash
    cp seq_folder/*cds.fa cds.fa
    cp seq_folder/*pro.fa pro.fa

    sed -i 's/\.protein_longest//g' species_order.txt
    cat species_order.txt | tr "\t " "\n" > so2.txt
    sp_number=$(cat so2.txt | wc -l)
    
    for s in *fa
    do
    cat \$s | grep ">" > patterns.txt
    for ((i=1; i<=sp_number; i++)); do
    species=$(awk -v row="${i}" 'NR==row' so2.txt)
    seq=$(awk -v row="${i}" 'NR==row' patterns.txt)
    echo $species
    echo $seq
    
    sed -i "s/$seq/$

    done

    done 
    


    prank -d=pro.fa -t=$params.tree -o=output_file -F -showxml



    """
}
