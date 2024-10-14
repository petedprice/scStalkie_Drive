nextflow.enable.dsl=2

include { BWA_index } from './modules/gatk/BWA_index.nf'
include { BWA } from './modules/gatk/BWA.nf'
include { MarkDuplicates } from './modules/gatk/MarkDuplicates.nf'
include { SplitNCigarReads } from './modules/gatk/SplitNCigarReads.nf'
include { pre_HaplotypeCaller } from './modules/gatk/pre_HaplotypeCaller.nf'
include { CollectWgsMetrics } from './modules/gatk/CollectWgsMetrics.nf'
include { HaplotypeCaller } from './modules/gatk/HaplotypeCaller.nf'
include { GenomicsDBImport } from './modules/gatk/GenomicsDBImport.nf'
include { GenotypeGVCF } from './modules/gatk/GenotypeGVCF.nf'
include { VariantFiltration } from './modules/gatk/VariantFiltration.nf'
include { MergeVcfs } from './modules/gatk/MergeVcfs.nf'
include { tabix } from './modules/tabix.nf'
include { vcftools_filt } from './modules/vcftools_filt.nf'
include { PopGenome_gff as PopGenome_gff } from './modules/PopGenome_gff.nf'

include { MKT_vcftools } from './modules/MKT_prep/MKT_vcftools.nf'
include { snpeff_eff } from './modules/MKT_prep/snpeff_eff.nf'


workflow {

	if (params.bwa_index == 'yes'){
	BWA_index()
	}

	samples_ch = Channel
		.fromPath(params.sample_info)
		.splitCsv()
		.map {row -> row[0]}

	contigs_ch = Channel
                .fromPath(params.chromsizes)
                .splitCsv()
                .map {row -> tuple(row[0],row[1], row[2])}

	aligned=BWA(samples_ch)


	cleaned=MarkDuplicates(aligned)
	splitNCRed=SplitNCigarReads(cleaned)
	preHC=pre_HaplotypeCaller(splitNCRed.combine(contigs_ch))
	VCed=HaplotypeCaller(preHC)

	All_VCFs=VCed
		.groupTuple( by:[1,2,3] )

	DBIed=GenomicsDBImport(All_VCFs)	
	GVCFed=GenotypeGVCF(DBIed)

	All_VCFs2=GVCFed
                .groupTuple( by:[1,2,3] )

	merged=MergeVcfs(All_VCFs2)

	filtered1=VariantFiltration(merged)
	filtered2=vcftools_filt(filtered1)


	VCF_tabixed=tabix(filtered2)	

	ctg_length_ch = Channel
                .fromPath(params.chromsizes)
                .splitCsv()
                .map {row -> tuple(row[0],row[3])}
		.unique()
		.view()

	//gff_PopGenome=PopGenome_gff(VCF_tabixed.combine(ctg_length_ch, by:0).view())


        annotated_synnon=snpeff_eff(filtered2)
	MKT_vcftooled=MKT_vcftools(annotated_synnon)

}




