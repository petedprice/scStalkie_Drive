process get_refs {
    cpus = 1
    memory = '4 GB'
    time = '4h'


    input:
    tuple val(species), val(ref)
    output:
    tuple val(species), file("${species}.gff"), file("${species}.protein.faa"), file("${species}.cds.fna")
    
    script:
    """
    #!/bin/bash
    echo $species
    echo $ref
    rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt ./
    grep -E ${ref} assembly_summary_refseq.txt | cut -f 20 > ftp_links.txt
    awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=\$0;asm=\$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_gff_files.sh
    awk 'BEGIN{FS=OFS="/";filesuffix="cds_from_genomic.fna.gz"}{ftpdir=\$0;asm=\$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_cds_files.sh
    awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz "}{ftpdir=\$0;asm=\$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_pro_files.sh

    source download_gff_files.sh
    source download_cds_files.sh
    source download_pro_files.sh

    gunzip *gz
    mv *genomic.gff ${species}.gff
    mv *cds_from_genomic.fna ${species}.cds.fna
    mv *protein.faa ${species}.protein.faa

   

    """
}
