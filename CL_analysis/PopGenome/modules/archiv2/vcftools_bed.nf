process vcftools_bed {

    maxRetries 4
    cpus = { 4 * task.attempt }
    memory = { 16.GB * task.attempt }
    time = {4.hour}

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
    
    if [ \$bedtype == Xgene ]; then
	echo "Removing genes from VCF"
    	vcftools --gzvcf ${contig}_filt.vcf.gz --exclude-bed $bed --recode --keep-INFO-all --stdout |sed 's/\\.:0,0:0/\\.|.:0,0:0/g'| bgzip -c > ${contig}_filt_type\$bedtype.vcf.gz
    elif [ \$bedtype == NCXG ]; then
        echo "Removing genes but keeping introns from VCF"
        vcftools --gzvcf ${contig}_filt.vcf.gz --exclude-bed $bed --recode --keep-INFO-all --stdout |sed 's/\\.:0,0:0/\\.|.:0,0:0/g'| bgzip -c > ${contig}_filt_type\$bedtype.vcf.gz
    elif [ \$bedtype == NC ]; then
	echo "Keeping only non coding regions"
	vcftools --gzvcf ${contig}_filt.vcf.gz --bed $bed --recode --keep-INFO-all --stdout |sed 's/\\.:0,0:0/\\.|.:0,0:0/g'| bgzip -c > ${contig}_filt_type\$bedtype.vcf.gz
    else 
	echo "doing something else, i think keeping a set of snps"
	vcftools --gzvcf ${contig}_filt.vcf.gz --positions keep_sites.txt --recode --keep-INFO-all --stdout |sed 's/\\.:0,0:0/\\.|.:0,0:0/g'| bgzip -c > ${contig}_filt_type\$bedtype.vcf.gz
    fi


    tabix -p vcf ${contig}_filt_type\$bedtype.vcf.gz

    """
}
