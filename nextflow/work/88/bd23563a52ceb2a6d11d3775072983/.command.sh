#!/bin/bash
/users/bop20pp/software/cellranger-7.2.0/bin/cellranger mkgtf 	/mnt/parscratch/users/bop20pp/Drive/ref_genome/OMAstalkie//Tdal_ST_mt.gtf 	cr_Tdal_ST_mt.gtf

/users/bop20pp/software/cellranger-7.2.0/bin/cellranger mkref 	--genome=stalkie_cellranger_reference 	--fasta=/mnt/parscratch/users/bop20pp/Drive/ref_genome/OMAstalkie//Tdal_ST_mt.fna 	--genes=cr_Tdal_ST_mt.gtf
