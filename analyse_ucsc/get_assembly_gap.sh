# this file just convert the raw assembly gap data into bed data so that we can run overlap enrichment on 
raw_assembly_gap_fn=/u/home/h/havu73/project-ernst/data/ucsc/mm10/gap.txt.gz
output_fn=/u/home/h/havu73/project-ernst/data/ucsc/mm10/clean_bed/clean_assembly_gap_mm10.bed
zcat $raw_assembly_gap_fn | awk -F'\t' 'BEGIN {OFS = "\t"}{print $2,$3,$4}' > $output_fn
gzip -f $output_fn