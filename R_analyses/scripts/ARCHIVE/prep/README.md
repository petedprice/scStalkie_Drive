## R Libraries
optparse, SCINA, Seurat, tidyverse, Matrix, scales, cowplot, RCurl, stringr, ggpubr, future, future.apply, DoubletFinder, remotes


## First create a conda environment with the relevant packages 
conda create -n seurat -c bioconda -c conda-forge r-seurat r-base=4.2.3 r-optparse r-tidyverse r-matrix r-scales r-cowplot r-rcurl r-stringr r-ggpubr r-scina r-remotes

Activate the conda env
```
source activate seurat
#or
conda activate seurat
```
open an R window and install doubletfinder using 

```{r}
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
```
Your envinroment should now be ready to go.


## 1.filtering.R

#### To run (example): 
```
source activate seurat

Rscript 1.filtering.R \
	-m /home/bop20pp/software/MeioticDrive2022/ \
	-d /fastdata/bop20pp/MeioticDrive2022/R/data/filtered \
	-o /fastdata/bop20pp/MeioticDrive2022/R/
```
#### Commands

-m Where the github respository is located \
-d Where your filtered cellranger data is stored \
-o Where you want your output data (plots and RData) to be saved

#### Description 
Filters 10X data using  
1. min.cells=3 (min number of cells supporting expression of a gene)
2. min.features=200 (Cells expressing at least 200 features)
3. mitoRatio<0.05 (Ratio of transcripts from mitochondrial origin)

#### Output
Various plots and RData titled filtered_seurat.RDatawhich contains the seurat object filtered_seurat

## 2.Integrate_normalise.R

#### To run (example): 
```
source activate seurat

Rscript 2.Integrate_normalise.R \
	-d /fastdata/bop20pp/MeioticDrive2022/R/outdata/filtered_seurat.RData \
	-o /fastdata/bop20pp/MeioticDrive2022/R/ \
	-t 1 \
	-s /home/bop20pp/software/MeioticDrive2022/R_analyses/data/samples.txt \
	-c /home/bop20pp/software/MeioticDrive2022/R_analyses/data/cell_cycle_markers_complete.csv
```
#### Commands
-d Location of RData containing filtered seurat object from (1) \
-o Where you want your output data (plots and RData) to be saved \
-t number of threads (default=1) \
-s txt file containing samples to keep in the analysis (Example see R_analyses/data/samples.txt) \
-c csv file containing cell cycle markers (Example see R_analyses/data/cell_cycle_markers_complete.csv) 

#### Description  
1. Removes doublets 
2. SCT normalises, regressing out cell cycle scores, MT ratio and nUMI. 
3. Integration
4. PCA, UMAP, TSNE, Cluster ID with varying resolutions (c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,  2.5, 3) using first 40 PCs (dims=1:40)

#### Output
Various plots and RData titled integrated_seurat.RData  which contains the seurat object seurat_integrated

## 3.cell_type_ID.R
#### To run (example): 
```
source activate seurat

Rscript 3.cell_type_ID.R \
	-d /fastdata/bop20pp/MeioticDrive2022/R/outdata/integrated_seurat.RData \
	-o /fastdata/bop20pp/MeioticDrive2022/R/ \
	-t 2 \
	-l /home/bop20pp/software/MeioticDrive2022/R_analyses/data/ortholog_table.txt 
```
#### Commands
-d Location of RData containing integrated, clustered seurat object from (2) \
-o Where you want your output data (plots and RData) to be saved \
-t number of threads (default=1) \
-l Location of txt file matrix containing ortholog information including marker genes \(Example see R_analyses/data/ortholog_table.txt) \
-s Which marker source to use (i.e. which column in ortholog_table.txt to use, refering to particular paper markers were sourced from) (default=comp_clusters, combination of Witt et al 2019 compiled clusters and the single-cell atlas)

#### Description  
1. Sets deafult resolution to 0.4 for clustering
2. Identifies clusters using markers from Witt et al 2019 extended list using SCINA and sctype (headers= customclassif and scina_labels) respectively. 

#### Output
Various plots and RData titled ${MARKER_SOURCE}_marker_seurat.RData  which contains the seurat object seurat_marker

