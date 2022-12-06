### This is a rather old version of the code. Because it defines the non coding region as being chrom length from ucsc subtract the cds regions. And this only produces data for the hg38 version.

genecode_raw_data_fn=/u/home/h/havu73/project-ernst/data/GENCODE/mouse/gencode.vM25.annotation.gtf.gz
genecode_coding_fn=/u/home/h/havu73/project-ernst/data/GENCODE/mouse/clean_bed/CDS_vM25.bed
# first, get a bed file of the coding region of the genome
num_comment_line=5 # number of lines that are comments about the data set and not really the tab separated data lines. Therefore, when we take out data, we should ignore these first num_comment_line lines
# In the original data file: 
# first column: chomosome name
# fourth column: genomic start location
# fifth column: genomic end location
# third column: feature type (we will filter out CDS: coding sequence feature type)
zcat $genecode_raw_data_fn | awk -v nIgnoreLine=$num_comment_line '{if (NR>nIgnoreLine) print $0}' | awk -F'\t' -v cds="CDS" '{if ($3==cds) print $1"\t"$4"\t"$5"\t"$3}' > $genecode_coding_fn
echo "Done getting bed file of CDS from genecode data"

# next, we rearrange the file of genecode_coding_fn so that genes of the same chromosome are clustered together and they are also arranged based on increasing coordinates
rearrange_bed_file_code=/u/home/h/havu73/project-ernst/source/chromHMM_utilities_unpublished/reorganize_bed_file.sh
rearranged_temp_cds_fn=/u/home/h/havu73/project-ernst/data/GENCODE/mouse/clean_bed/CDS_vM25_rearranged.bed
$rearrange_bed_file_code $genecode_coding_fn $rearranged_temp_cds_fn
echo "Done rearranging the bed file of CDS region in the genome"
echo ""
rearranged_temp_cds_fn=${rearranged_temp_cds_fn}.gz # change to .gz file after running the code
# next, we get rid of overlapping bases of the coding region
eradicate_overlapping_bases_code=/u/home/h/havu73/project-ernst/source/chromHMM_utilities_unpublished/eradicate_overlaping_bases_bed.py
non_overlap_cds_fn=/u/home/h/havu73/project-ernst/data/GENCODE/mouse/clean_bed/CDS_vM25_no_overlap.bed.gz
python $eradicate_overlapping_bases_code $rearranged_temp_cds_fn $non_overlap_cds_fn
echo "Done getting rid of overlapping bases of the cds bed file"
echo ""
# next, clean up files
# gzip $non_overlap_cds_fn 

