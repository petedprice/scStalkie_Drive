process longest_isoform_agat {

    cpus = 4
    memory = '128 GB'
    time = '4h'

    label 'agat'

    input:
    tuple val(species), val(USE) //, file("${species}.gff"), file("${species}.protein.faa"), file("${species}.cds.fna")
 
    output:
    tuple file("${species}.cds_longest.fna"), file("${species}.protein_longest.faa")

    script:
    """
    #!/bin/bash

    agat_sp_keep_longest_isoform.pl -gff ${launchDir}/refs/${species}.gff -o ${species}_longest.gff
    agat_sp_extract_sequences.pl -g ${species}_longest.gff -f ${launchDir}/refs/${species}.fna -t cds -p -o ${species}.protein_longest.faa
    agat_sp_extract_sequences.pl -g ${species}_longest.gff -f ${launchDir}/refs/${species}.fna -t cds -o ${species}.cds_longest_tmp.fna

    awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < ${species}.cds_longest_tmp.fna | tail -n +2 > ${species}.cds_longest.fna


    """
}
