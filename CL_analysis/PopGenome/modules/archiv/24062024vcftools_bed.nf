process vcftools_bed {

    maxRetries 4
    cpus = { 16 * task.attempt }
    memory = { 256.GB * task.attempt }
    time = {12.hour}

    label 'vcftools'

    tag {'vcftools_bed_' + contig }

    input:
    tuple val(vcf), val(contig), val(ctg_len), file(bed), file("${contig}_filt.vcf.gz")

    output:
    tuple val(contig), val(ctg_len), env(bedtype), file('*type*vcf.gz'), file('*type*vcf.gz.tbi')    
    
    script:
    """
    #!/bin/bash
    
    basename=\$(basename $bed .bed)
    bedtype="\${basename##*_}"

    echo \$bedtype
    cat $bed |cut -f 1,3 >> keep_sites.txt
    
    if [ \$bedtype==xGene ]; then
    	vcftools --gzvcf ${contig}_filt.vcf.gz --bed $bed --recode --keep-INFO-all --stdout | bgzip -c > ${contig}_filt_type\$bedtype.vcf.gz
    else
    	#zcat ${contig}_filt.vcf.gz | grep -F -f keep_sites.txt | bgzip -c > ${contig}_filt_type\$bedtype.vcf.gz 
	vcftools --gzvcf ${contig}_filt.vcf.gz --positions keep_sites.txt --recode --keep-INFO-all --stdout | bgzip -c > ${contig}_filt_type\$bedtype.vcf.gz
    fi

    tabix -p vcf ${contig}_filt_type\$bedtype.vcf.gz

    """
}
