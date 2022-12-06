# this file just convert the raw rmsk data into bed data so that we can run overlap enrichment on 
raw_rmsk_fn=/u/home/h/havu73/project-ernst/data/ucsc/mm10/rmsk.txt.gz
output_fn=/u/home/h/havu73/project-ernst/data/ucsc/mm10/clean_bed/repeats.bed
zcat $raw_rmsk_fn | awk -F'\t' 'BEGIN {OFS = "\t"}{print $6,$7,$8}' > $output_fn
gzip -f $output_fn