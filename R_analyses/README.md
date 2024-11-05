## Rscripts to run analysis downstream analysis of output from nextflow pipeline 

Initial input includes a Seurat object for all samples combined. 

## Scripts ordering


### Orthologs.R

Create consensus set of orthologs between those identified by UCL and results from running orthofinder against the drosophila reference genome. 


### 1_Joint_Celltypes_Visualise_data_full.R

Script to run initial cell type identification and filtering of low confidence cells and clusters 


### 2_Joint_Celltypes_Visualise_data_subset.R

Script to re-run cell type identification on the subsetted data
Outputs Seurat object used in all other scripts

### DGE_DC_DA.R

Differential gene expression and dosage compensation analysis

### Differential_Abundance.R

Differential abundances between ST and SR samples across cell types 


### Ploidy_check.R

Check ploidy of cells for supplementary results 

### Trajectory_analaysis_cluster.R

Initial run of trajectory analysis for running on the cluster to find the correct number of knots for the GAM 


### Trajectory_Analysis.R

Trajectory analysis on object from Trajectory_analaysis_cluster.R


### SUP_TABLES.R

Create supplementary tables

### Plots.R

Plots for manuscript

### MSL_plots.R

Plotting of the expression of MSL genes 


