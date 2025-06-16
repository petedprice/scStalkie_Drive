nextflow.enable.dsl=2

include { paml_tree } from './modules/paml_tree.nf'
include { paml_branch_tree } from './modules/paml_branch_tree.nf'

//include { swamp_branches } from './modules/swamp_branches.nf'


include { get_refs } from './modules/get_refs.nf'
include { longest_isoform_agat } from './modules/longest_isoform_agat.nf'

include { orthofinder } from './modules/orthofinder.nf'
include { ortho_cds } from './modules/ortho_cds.nf'

include { prank_allign } from './modules/prank_allign.nf'
include { format_fastas } from './modules/format_fastas.nf'
include { prank_phy } from './modules/prank_phy.nf'

include { remove_gaps } from './modules/remove_gaps.nf'
include { remove_Ns } from './modules/remove_Ns.nf'

include { mod0_paml } from './modules/mod0_paml.nf'
include { m1avsm2a } from './modules/m1avsm2a.nf'
include { branch_model } from './modules/branch_model.nf'
include { comp_paml_models } from './modules/comp_paml_models.nf'

include { swamp_test } from './modules/swamp_test.nf'
include { swamp_final } from './modules/swamp_final.nf'

include { get_omega } from './modules/get_omega.nf'


workflow {
    //Channels species name and reference name
    if (params.swamp == "checks"){
	    swamp_test_parameters=Channel
        	.fromPath(params.swamp_test_params)
        	.splitCsv()
        	.map {row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])}
        	.unique()
		.view()
    } else if(params.swamp == "yes"){
	    swamp_final_parameters=Channel
                .fromPath(params.swamp_final_params)
                .splitCsv()
                .map {row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])}
                .unique()
		.view()

    }


    tree=Channel.fromPath(params.tree)
		.view()

   for( line : file(params.tree).readLines() ) {
	println line
    }


    branch_tree=Channel
	        .fromPath(params.branch_tree, type: 'file')
		.view()

    trees=tree.combine(branch_tree)
		.view()

    pamled_tree=paml_tree(tree)

    pamled_branch_tree=paml_branch_tree(trees)


    species_ch=Channel
		.fromPath(params.metadata)
		.splitCsv()
		.map {row -> tuple(row[0], row[1])}
		.unique()
		.view()

    longest=longest_isoform_agat(species_ch).collect()
    orthofound=orthofinder(longest)
    of_cds=ortho_cds(orthofound)
	.transpose()
	
    species_formatted=format_fastas(of_cds)

    pranked_alligned=prank_allign(species_formatted.combine(tree))
    No_gaps=remove_gaps(pranked_alligned)
    phyed=prank_phy(No_gaps)

    

    if (params.swamp !="no"){
	    mod0ed=mod0_paml(phyed.combine(pamled_tree))
    }


    if(params.swamp == "checks"){
	   swamped=swamp_test(mod0ed.combine(swamp_test_parameters))
		.transpose()
	   M1aM2aed=m1avsm2a(swamped.combine(pamled_tree))
                .groupTuple(by: 1)
		.unique()

	   top_mod=comp_paml_models(M1aM2aed)
           //top_algns=return_top_algns(top_mod.combine(swamped.groupTuple(by: 0)))

    } else if(params.swamp == "yes"){
	   swamped=swamp_final(mod0ed.combine(swamp_final_parameters))
	   NoNs=remove_Ns(swamped)
           //M1aM2aed=m1avsm2a(NoNs.combine(pamled_tree))
	//	.groupTuple(by: 1)
           //top_mod=comp_paml_models(M1aM2aed)
           

           branch_modeled=branch_model(NoNs.combine(pamled_branch_tree))
                .groupTuple(by: 1)
	   //branch_omegas=get_omega(branch_modeled)

/*	   
	   combined=top_mod
		.combine(swamped
			.groupTuple(by: 0))
		.combine(NoNs
			.groupTuple(by: 0))
		.view()
 	   
*/
//top_algns=return_top_algns(top_mod.combine(swamped.groupTuple(by: 0)).combine(NoNs.groupTuple(by: 0)

    } else if(params.swamp == "no"){
           M1aM2aed=m1avsm2a(phyed.combine(pamled_tree))
                .groupTuple(by: 1)
           top_mod=comp_paml_models(M1aM2aed)	     

    }

}




