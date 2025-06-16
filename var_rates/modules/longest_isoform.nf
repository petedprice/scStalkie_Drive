process longest_isoform {

    cpus = 4
    memory = '128 GB'
    time = '4h'

    //label 'python'

    input:
    tuple val(species), val(USE) //, file("${species}.gff"), file("${species}.protein.faa"), file("${species}.cds.fna")
 
    output:
    tuple file("${species}.cds_longest.fna"), file("${species}.protein_longest.faa")

    script:
    """
    #!/bin/bash

    cp ${launchDir}/refs/${species}.protein.faa .
    python ${baseDir}/software/OrthoFinder_source/tools/primary_transcript.py ${species}.protein.faa
    mv primary_transcripts/${species}.protein.faa ${species}.protein_longest.faa

    cat ${species}.protein_longest.faa | grep '>' | cut -f1 -d " " | cut -c 2- > longest_proteins.txt
    awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < ${launchDir}/refs/${species}.cds.fna | tail -n +2 > singleline_cds.fasta
    grep -F -A 1 --no-group-separator -f longest_proteins.txt singleline_cds.fasta > ${species}.cds_longest.fna



    """
}
