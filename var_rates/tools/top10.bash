mkdir top_genes
cat $1 | head -11 | tail -n +2 | cut -f5 -d ' ' > top_genes.txt
for i in $(cat top_genes.txt) 
do 
echo $i
cp unmasked_allignments/${i}* top_genes
cp swamp_masked_allignments/${i}* top_genes
for p in top_genes/*${i}*
do 
Rscript ~/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/scripts/phylip_to_fa.R $p
done
done
