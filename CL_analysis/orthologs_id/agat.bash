source  activate agat

agat_sp_keep_longest_isoform.pl -gff $1/$2.gff -o $1/${2}_longest.gff

agat_sp_extract_sequences.pl -gff $1/${2}_longest.gff --fasta $1/$2.fasta -t cds -o $1/${2}_longest.fasta

