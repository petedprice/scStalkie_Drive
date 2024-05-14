#!/bin/bash

#$ -l h_rt=8:0:0

#$ -l rmem=32G

#$ -pe smp 8

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/scRNAseq/wdir

source /usr/local/extras/Genomics/.bashrc

seq1="$(basename $1)"
seq2="$(basename $2)"
seq1=${seq1%.fasta}
seq2=${seq2%.fasta}

dir1="$(dirname $1)"
dir2="$(dirname $2)"

#makeblastdb -in $1 -input_type fasta -dbtype nucl -out ${seq1} -parse_seqids
#makeblastdb -in $2 -input_type fasta -dbtype nucl -out ${seq2} -parse_seqids

makeblastdb -in $1 -input_type fasta -dbtype nucl -parse_seqids
makeblastdb -in $2 -input_type fasta -dbtype nucl -parse_seqids

echo "running blast one"

blastn -evalue 1e-10 \
	-task blastn \
	-db $1 \
	-query $2 \
	-out ${dir1}/${seq1}.${seq2}.bla \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'

echo "running blast two"


blastn -evalue 1e-10 \
        -task blastn \
        -db $2 \
        -query $1 \
        -out ${dir1}/${seq2}.${seq1}.bla \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'

