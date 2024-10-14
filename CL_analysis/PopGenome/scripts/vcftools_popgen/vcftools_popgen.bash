vcftools --gzvcf $1 --weir-fst-pop $2 --weir-fst-pop $3 --fst-window-size 100000 --fst-window-step 25000 --out fst 
vcftools --gzvcf $1 --window-pi 100000 --window-pi-step 25000 --keep $2 --out ST_pi
vcftools --gzvcf $1 --window-pi 100000 --window-pi-step	25000 --keep $3	--out SR_pi
vcftools --gzvcf $1 --TajimaD --keep $2 --out ST_taj
vcftools --gzvcf $1 --TajimaD --keep $3 --out SR_taj
