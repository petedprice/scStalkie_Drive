gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /mnt/parscratch/users/bop20pp/Drive/ref_genome/Tdal_ST_mt.fna \
   -V gendb://${1}_drive_genomicsdb \
   -O /mnt/parscratch/users/bop20pp/Drive/PopGenome/gatk/${1}_ST_SR.vcf.gz
