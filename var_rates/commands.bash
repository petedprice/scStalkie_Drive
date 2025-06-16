module load Nextflow/22.04.0
module load Anaconda3/2022.05

nextflow run /users/bi1pp/personal_git/scStalkie_Drive/var_rates/main.nf \
	--metadata /users/bi1pp/personal_git/scStalkie_Drive/var_rates/metadata.csv \
	-resume



