# Pop genetics and variant calling pipeline etc

User must have nextflow installed on their system to run. 
To run, the below example command may be used

```
module load Nextflow/22.04.0

nextflow run /users/bop20pp/personal_git/scStalkie_Drive/CL_analysis/PopGenome/main.nf \
	--fasta /mnt/parscratch/users/bop20pp/Drive/ref_genome/Tdal_ST_mt.fna \
        --gff /mnt/parscratch/users/bop20pp/Drive/ref_genome/Tdal_ST_mt.gff \
        --sample_info  /users/bop20pp/personal_git/scStalkie_Drive/CL_analysis/PopGenome/sample_info.csv \
        --reads /mnt/parscratch/users/bop20pp/Drive/PopGenome/data/sequencing_data/StalkEyesFlies_Sep2022/Trimmed/Sample_ \
	--chromsizes /users/bop20pp/personal_git/scStalkie_Drive/CL_analysis/PopGenome/chromsizes_bin.csv \
	--snpeff_db /mnt/parscratch/users/bop20pp/Drive/PopGenome/data/snpeff_data/stalkie \
	-resume

```


main.nf is the main pipeline script 


it calls modules from the folder 'modules' or 'modules/gatk'



some of these modules call scripts from the 'scripts' folder 


you may need to modify the nextflow.config file to change directories for where containers are stored etc



enjoy!
