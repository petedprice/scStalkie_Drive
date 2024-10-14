nextflow.enable.dsl=2

include { tabix } from './modules/tabix.nf'
include { vcftools_filt } from './modules/vcftools_filt.nf'
include { PopGenome_vcf as PopGenome_gff } from './modules/PopGenome_gff.nf'

workflow {





	contig_ch=Channel.fromPath(params.metadata)
		.splitCsv()
		.map {row -> tuple(row[0], row[1], row[2])}
		.unique()

	VCF_filtered=vcftools_filt(contig_ch)	
	VCF_nobeded=tabix(VCF_filtered)	
	gff_PopGenome=PopGenome_gff(VCF_nobeded)
 
}




