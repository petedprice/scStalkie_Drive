module load Nextflow/22.04.0


nextflow run /users/bi1pp/personal_git/scStalkie_Drive/nextflow/main.nf \
	--fasta_dir /mnt/parscratch/users/bi1pp/Drive/ref_genome \
	--gtf_dir /mnt/parscratch/users/bi1pp/Drive/ref_genome \
	--cellranger /users/bi1pp/software/cellranger-7.2.0/bin/cellranger \
	--metadata /users/bi1pp/personal_git/scStalkie_Drive/nextflow/metadata_mt.csv \
	--read_dir /mnt/parscratch/users/bi1pp/Drive/reads \
	--cellcycle_markers /users/bi1pp/personal_git/scStalkie_Drive/nextflow/data/cell_cycle_markers.csv \
	--celltype_markers /users/bi1pp/personal_git/scStalkie_Drive/nextflow/data/markers_elife_plus_extras.csv \
	--mt_genes /users/bi1pp/personal_git/scStalkie_Drive/nextflow/data/mt_genes.txt \
	-resume




