# Nextflow pipeline for running cell ranger and Seurat on 10X sequencing data 

## Commands 
module load Nextflow/22.04.0

```
nextflow run /users/bop20pp/personal_git/scStalkie_Drive/nextflow/main.nf \
	--fasta_dir /mnt/parscratch/users/bop20pp/Drive/ref_genome \
	--gtf_dir /mnt/parscratch/users/bop20pp/Drive/ref_genome \
	--cellranger /users/bop20pp/software/cellranger-7.2.0/bin/cellranger \
	--metadata /users/bop20pp/personal_git/scStalkie_Drive/nextflow/metadata_mt.csv \
	--read_dir /mnt/parscratch/users/bop20pp/Drive/reads \
	--cellcycle_markers /users/bop20pp/personal_git/scStalkie_Drive/nextflow/data/cell_cycle_markers.csv \
	--celltype_markers /users/bop20pp/personal_git/scStalkie_Drive/nextflow/data/markers_elife_plus_extras.csv \
	-resume \
	-with-dag flowchat.png \
	-with-trace
```

fasta_dir: directory where fasta reference is stored 

gtf_dir: directory where gtf reference is stored 

cellranger: path to cellranger 

metadata: sample metadata 

read_dir: path to reads

cellcycle_markers: markers of cell cycle 

celltype_markers: celltype markers (will give preliminary cell types but these are not used in final manuscript) 








# Notes 

Singularity/apptainer tips
Make sure in your bashrc to ensure enough space to pull containers. 
export SINGULARITY_CACHEDIR=/fastdata/bop20pp/apptainer
export NXF_SINGULARITY_CACHEDIR=/fastdata/bop20pp/apptainer
export APPTAINER_CACHEDIR=/fastdata/bop20pp/apptainer
export NXF_APPTAINER_CACHEDIR=/fastdata/bop20pp/apptainer
export SINGULARITY_TMPDIR=/fastdata/bop20pp/apptainer
export APPTAINER_TMPDIR=/fastdata/bop20pp/apptainer

