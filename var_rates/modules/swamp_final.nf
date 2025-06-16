process swamp_final {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'swamp_final' + '_' + og }


    publishDir 'swamp_masked_allignments_final', mode: 'copy', overwrite: true, pattern: '*.phy'

    label 'python'

    input:
    tuple file("${og}_swamp_analysis"), val(og), file("${og}_SWAMP_BRANCHES.txt"), val(w1), val(t1), val(w2), val(t2), val(w3), val(t3), val(w4), val(t4) 

    output:
    tuple file('*.phy'), val(og)
    
    script:
    """
    #!/bin/bash
    mkdir tmp_swamp
    mkdir final_phy
    cp ${og}_swamp_analysis/* tmp_swamp


    if [ $w1 != "null" ]; then 
    python2.7 ${baseDir}/software/SWAMP/SWAMP.py -i tmp_swamp -b ${og}_SWAMP_BRANCHES.txt \
	-t $t1 -w $w1 -m 100 > ${og}_swamp_t${t1}w${w1}.txt
    mv tmp_swamp/${og}_codon.nogaps_masked.phy tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}.phy
    cp tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}.phy final_phy
    rm tmp_swamp/${og}_codon.nogaps.phy
    fi

    if [ $w2 != "null" ]; then
    python2.7 ${baseDir}/software/SWAMP/SWAMP.py -i tmp_swamp -b ${og}_SWAMP_BRANCHES.txt \
        -t $t2 -w $w2 -m 100 > ${og}_swamp_t${t1}w${w1}_t${t2}w${w2}.txt
    rm final_phy/*
    mv tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_masked.phy tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}.phy
    rm tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}.phy
    cp tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}.phy final_phy/
    fi   

    if [ $w3 != "null" ]; then
    python2.7 ${baseDir}/software/SWAMP/SWAMP.py -i tmp_swamp -b ${og}_SWAMP_BRANCHES.txt \
        -t $t3 -w $w3 -m 100 > ${og}_swamp_t${t1}w${w1}_t${t2}w${w2}_t${t3}w${w3}.txt
    rm final_phy/*
    mv tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}_masked.phy tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}_t${t3}w${w3}.phy
    cp tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}_t${t3}w${w3}.phy final_phy   
    rm tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}.phy
    fi

    if [ $w4 != "null" ]; then
    python2.7 ${baseDir}/software/SWAMP/SWAMP.py -i tmp_swamp -b ${og}_SWAMP_BRANCHES.txt \
	-t $t3 -w $w3 -m 100 > ${og}_swamp_t${t1}w${w1}_t${t2}w${w2}_t${t3}w${w3}_t${t4}w${w4}.txt
    rm final_phy/*
    mv tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}_t${t3}w${w3}_masked.phy tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}_t${t3}w${w3}_t${t4}w${w4}.phy
    cp tmp_swamp/${og}_codon.nogaps_t${t1}w${w1}_t${t2}w${w2}_t${t3}w${w3}_t${t4}w${w4}.phy final_phy
    fi

    cp final_phy/*phy .        

    """
}

