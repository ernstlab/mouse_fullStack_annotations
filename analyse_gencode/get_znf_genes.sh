#raw_genecode_data_fn=/u/home/h/havu73/project-ernst/data/GENCODE/gencode.v30lift37.annotation.gtf.gz
raw_genecode_data_fn=/u/home/h/havu73/project-ernst/data/GENCODE/mouse/gencode.vM25.annotation.gtf.gz
# first get any lines with zfp in the lines and then print columns 1, 4, 5, 9 because those correspond to chromosome, start, end, gene annotations
step1_output_fn=${PWD}/zfp_step1
zcat $raw_genecode_data_fn | grep "Zfp" | awk -F'\t' 'BEGIN {OFS="\t"}{print $1,$4,$5,$9}' > $step1_output_fn # I tried some variations of zfp such as zfp, zfp, zfp, and the zfp is the one that worked for this particular dataset

# second rearrange the genes so that they are all ordered along the genome
rearrange_bed_file_code=/u/home/h/havu73/project-ernst/source/chromHMM_utilities_unpublished/reorganize_bed_file.sh
step2_output_fn=${PWD}/zfp_step2
$rearrange_bed_file_code $step1_output_fn $step2_output_fn
step2_output_fn=${step2_output_fn}.gz

# lastly remove any overlapping bases between the zfp genes' annotations
eradicate_overlapping_bases_bed_code=/u/home/h/havu73/project-ernst/source/chromHMM_utilities_unpublished/eradicate_overlaping_bases_bed.py
# CHANGE THIS IF YOU WANT TO CHANGE BETWEEN HG19 TO HG38
#clean_zfp_bed_fn=/u/home/h/havu73/project-ernst/data/GENCODE/zfp_genes/hg19/zfp_genes_no_overlap.gz
clean_zfp_bed_fn=/u/home/h/havu73/project-ernst/data/GENCODE/mouse/clean_bed/zfp_genes_no_overlap.gz
python $eradicate_overlapping_bases_bed_code $step2_output_fn $clean_zfp_bed_fn

#rm $step1_output_fn
#rm $step2_output_fn
