module load Nextflow/22.04.0


nextflow run /users/bop20pp/personal_git/scStalkie_Drive/nextflow/main.nf \
	--fasta_dir /mnt/parscratch/users/bop20pp/Drive/ref_genome \
	--gtf_dir /mnt/parscratch/users/bop20pp/Drive/ref_genome \
	--cellranger /users/bop20pp/software/cellranger-7.2.0/bin/cellranger \
	--metadata /users/bop20pp/personal_git/scStalkie_Drive/nextflow/metadata_mt.csv \
	--read_dir /mnt/parscratch/users/bop20pp/Drive/reads \
	--cellcycle_markers /users/bop20pp/personal_git/scStalkie_Drive/nextflow/data/cell_cycle_markers.csv \
	--celltype_markers /users/bop20pp/personal_git/scStalkie_Drive/nextflow/data/markers_elife_plus_extras.csv \
	--mt_genes /users/bop20pp/personal_git/scStalkie_Drive/nextflow/data/mt_genes.txt \
	-resume \
	-with-trace




