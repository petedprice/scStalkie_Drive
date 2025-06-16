# Detecting Coding Sequence Evolution

A nextlow pipeline utilising PRANK, SWAMP, and PAML to detecting evolution on the coding sequence across a number of optional species

The project is oranigsed into [modules](https://github.com/petedprice/Avian_scRNAseq/tree/main/CL_analyses/nextflow/var_rates/modules), each one running a specific process/script. 
These are controlled in the [main.nf](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/main.nf) file. 
User specific controls can be altered in the [nextflow.config](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/nextflow.config) file. 
Likewise some [data](https://github.com/petedprice/Avian_scRNAseq/tree/main/CL_analyses/nextflow/var_rates/data) input is saved in the data files of which there are defaults that can be changed for user preferences.


User must have nextflow installed on their system to run. 
To run, the below example command may be used

```
module load Nextflow/22.04.0
module load Anaconda3/2022.10

nextflow run /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/main.nf \
	--metadata /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/metadata_full.csv \
	--tree /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/data/tree.txt \
	-resume \
	-with-dag flowchat.png \
	-with-trace
```

# INPUT 

**metadata**
- two column csv file with species name (seperated parts by underscore) in first column and NCBI RefSeq assembly in second column

**tree**
- This is your nexus tree for running orthofinder 
- This will be converted into a unrooted tree for PAML and branch file will be made from it for SWAMP 


**swamp_test_params**
- Parameters to run through when testing SWAMP masking 


**swamp_final_params**
- Final parameters to use with SWAMP. Masked regions will then be removed from the allignment 

**swamp**
- Either ```yes, no or checks```
- yes: will run SWAMP once with a single set of paremeters and output masked allignements with Ns and once Ns are removed. PAML will be run on the N removed allignmenets
- checks: will run SWAMP on all parameters combinations and output masked allignments. PAML will be run on allignements with Ns included 
- no: will not mask


# PIPELINE

## Core analaysis

 **[get_refs](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/get_refs.nf)**
- Script pulls down user defined CDS and proteome from list of species from NCBI
- To do:
		- Include alternate mode if CDS/proteome aren't available, i.e. predict ORF and use gene seq to generate orthologs)


**[longest_isoform](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/longest_isoform.nf)**
- Get longest isoform for each protein for each species 

**[orthofinder](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/orthofinder.nf)** 
- Creates one-to-one ortholog list for set of species used and returns a one-to-one proteome and transcriptome
- More info at the [Orthofinder Github](https://github.com/davidemms/OrthoFinder)

**[ortho_cds](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/ortho_cds.nf)**
- From the one-to-one transcriptome, generates an individual fasta for each orthogroup

**[format_fastas](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/format_fastas.nf)**
- For each orthogroup fasta for both protein and cds, rename gene names with species name 

**[prank_allign](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/prank_allign.nf)**
- For each orthogroup, allign CDS fasta
- For more information on [PRANK](http://wasabiapp.org/software/prank/)

**[remove_gaps](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/remove_gaps.nf)**
- Remove gaps in the allignment using [remove_gaps.R](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/scripts/remove_gaps.R)

**[prank_phy](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/prank_phy.nf)**
- Convert allignment into Phylip format for input into PAML

## Option 1: No masking (swamp == "no")

**[m1avsm2a](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/m1avsm2a.nf)**
- Run the M1A and M2A model in PAML from your unmasked allignment from prank

**[comp_paml_models](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/comp_paml_models.nf)**
- Using likliehood ratios compare fit of both models and output summary in single output for all genes 

## Option 2: check masking parameters (swamp == "checks")

**[mod0_paml](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/mod0_paml.nf)**
- Run model 0 in PAML, as the required input for SWAMP
- To do, add more information on this model

**[swamp_test](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/swamp_test.nf)**
- Run SWAMP on set of parameters selected by user in unheaded CSV file. See [example](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/data/SWAMP_PARAMS/SWAMP_TEST.txt)
- User can have up to 4 masking combinations, e.g. window size of 15 and number of non-synonymous differences as 2 as initial mask. Then layered on top a w of 10 and t of 7. Then onto a layer of.... and so on. 
- User can state as many parameter mixes as needed on new lines in the input file
- Fill blanks with null

**[m1avsm2a] (https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/m1avsm2a.nf)** 
- Site model comparisons: m1a is nearly neutral and m2a is positive selection
- for every combination and subset of masking, models will be run and compared
- Can alter the stock [file](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/data/PAML_CTLs/mod1a2a.ctl) to choose the models you want to compare. e.g. if you also want to run model 7 and 8 however model comparisons will vary downstream
- To do: include extra model parameter flexilbity 

**[comp_paml_models](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/comp_paml_models.nf)**

## Option 3: Final SWAMP masking (swamp == "yes")

**[swamp_final](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/swamp_final.nf)**
- Run SWAMP on final parameter set selected by user in unheaded CSV file. See [example](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/data/SWAMP_PARAMS/SWAMP_FINAL.txt)
- User can have up to 4 masking combinations, e.g. window size of 15 and number of non-synonymous differences as 2 as initial mask. Then layered on top a w of 10 and t of 7. Then onto a layer of.... and so on. 
- User only provide a single line
- fill blanks with null
- output to next process will only be complete masking, not the intermediate masking values. e.g. if you have two masking rounds, only the cumalative masking after the second round will be taken to the next step

**[remove_Ns](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/remove_Ns.nf)**
- Remove Ns from the allignment 


**[m1avsm2a](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/m1avsm2a.nf)**
- For masked allignment run m1a and m2a models

**[comp_paml_models](https://github.com/petedprice/Avian_scRNAseq/blob/main/CL_analyses/nextflow/var_rates/modules/comp_paml_models.nf)**


# TO DO 
- Push changes 
- Include branch model 
- Include M0 model on masked data
- Compile M0, M1aM2a data into the comp models script 


