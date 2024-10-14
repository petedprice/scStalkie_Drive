process MarkDuplicates {

    label 'gatk'
    errorStrategy 'retry'

    tag {'markduplicates_' + sample }

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    input:
    tuple val(sample), file("${sample}.sam")

    output:
    tuple val(sample), file("${sample}_marked_duplicates.bam")

    script:
    """
    #!/bin/bash
    
   gatk --java-options "-Xmx4g" SortSam \
	-I ${sample}.sam \
	-O ${sample}_sorted.sam \
	-SORT_ORDER coordinate

    gatk --java-options "-Xmx4g" MarkDuplicates \
	-I ${sample}_sorted.sam \
	-O ${sample}_marked_duplicates.bam \
	-M ${sample}_dup_metrics.txt \
	-ASSUME_SORTED true \
	--VALIDATION_STRINGENCY SILENT

    """
}

