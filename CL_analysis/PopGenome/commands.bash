module load Nextflow/22.04.0

nextflow run /users/bop20pp/personal_git/scStalkie_Drive/CL_analysis/PopGenome/main.nf \
	--fasta /mnt/parscratch/users/bop20pp/Drive/ref_genome/Tdal_ST_mt.fna \
        --gff /mnt/parscratch/users/bop20pp/Drive/ref_genome/Tdal_ST_mt.gff \
        --sample_info  /users/bop20pp/personal_git/scStalkie_Drive/CL_analysis/PopGenome/sample_info.csv \
        --reads /mnt/parscratch/users/bop20pp/Drive/PopGenome/data/sequencing_data/StalkEyesFlies_Sep2022/Trimmed/Sample_ \
	--chromsizes /users/bop20pp/personal_git/scStalkie_Drive/CL_analysis/PopGenome/chromsizes_bin.csv \
	--snpeff_db /mnt/parscratch/users/bop20pp/Drive/PopGenome/data/snpeff_data/stalkie \
	-resume

	#-resume


	#--metadata /users/bop20pp/personal_git/scStalkie_Drive/CL_analysis/PopGenome/contig_metadata.csv \
	#--gff /mnt/parscratch/users/bop20pp/Drive/ref_genome/Tdal_ST_mt.gtf \
	#--vcfs /mnt/parscratch/users/bop20pp/Drive/PopGenome/vcfs \
