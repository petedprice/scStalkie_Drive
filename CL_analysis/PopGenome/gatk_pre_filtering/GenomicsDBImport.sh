gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
       --sample-name-map /users/bop20pp/personal_git/scStalkie_Drive/CL_analysis/PopGenome/gatk_pre_filtering/cohort.sample_map \
      --genomicsdb-workspace-path /mnt/parscratch/users/bop20pp/Drive/PopGenome/gatk/${2}_drive_genomicsdb \
      --tmp-dir /mnt/parscratch/users/bop20pp/Drive/PopGenome/gatk/tmp \
	-L $1
