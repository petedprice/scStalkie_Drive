conda activate bedops 

input_file="$1"
output_file="${input_file%.vcf}.bed"

vcf2bed < "$input_file" | cut -f1,3,6,7 > tmp.vcf
cat tmp.vcf | grep -v "*" > tmp2.vcf
expand -t 1 tmp2.vcf > "$output_file"
