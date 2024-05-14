source activate seurat
Rscript 2.Integrate_normalise.R -d /fastdata/bop20pp/MeioticDrive2022/R/outdata/filtered_seurat.RData -o /fastdata/bop20pp/MeioticDrive2022/R/ -t 1 -s /home/bop20pp/software/MeioticDrive2022/R_analyses/data/samples.txt \
	-c /home/bop20pp/software/MeioticDrive2022/R_analyses/data/cell_cycle_markers_complete.csv
