source activate seurat
Rscript 3.cell_type_ID.R -d /fastdata/bop20pp/MeioticDrive2022/R/outdata/integrated_seurat.RData \
	-o /fastdata/bop20pp/MeioticDrive2022/R/ -t 2 \
	-l /home/bop20pp/software/MeioticDrive2022/R_analyses/data/ortholog_table.txt 	
