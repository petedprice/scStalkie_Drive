#!/bin/bash

mkdir contigs
cat /mnt/parscratch/users/bop20pp/Drive/ref_genome/OMAstalkie//Tdal_ST_mt.gtf | cut -f 1 | uniq > contigs.txt
for contig in $(cat contigs.txt)
do
echo $contig > ${contig}_stalkie.txt
done
rm contigs.txt
