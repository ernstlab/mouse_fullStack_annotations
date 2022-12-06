genecode_data_fn=/u/home/h/havu73/project-ernst/data/GENCODE/mouse/gencode.vM25.annotation.gtf.gz
pseudogene_fn=/u/home/h/havu73/project-ernst/data/GENCODE/mouse/clean_bed/gene_type_pseudogene_vM25_mm10.bed
#genecode_data_fn=/u/home/h/havu73/project-ernst/data/GENCODE/gencode.v30lift37.annotation.gtf.gz
#pseudogene_fn=/u/home/h/havu73/project-ernst/data/GENCODE/gene_type_pseudogene_v30lift37.bed

rm -f $pseudogene_fn
# echo -e "chromosome\tstart_bp\tend_bp\tfeature_type\tgene_type\ttranscript_type" > $pseudogene_fn
num_comment_line=5 # number of lines that are comments about the data set and not really the tab separated data lines. Therefore, when we take out data, we should ignore these first num_comment_line lines
# this line of code is old and it does not work for both hg19 and hg38 data. So, I replaced it to selecting genes (lines) that contains the pseudogene word in their annotaiton (grep pseudogene) zcat $genecode_data_fn | awk -v nIgnoreLine=$num_comment_line '{if (NR>nIgnoreLine) print $0}' | awk -F'\t' -v cds="CDS" '{if ($3==cds) print $0}' | awk -F'\t' '{gsub("; ","\t",$9); print $0}' | awk -F'[\t ]' 'BEGIN { OFS="\t" } {gsub(/"/,"",$0); if ($12 ~ /pseudogene/) print $0}' | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$11,$12}' >> $pseudogene_fn
zcat $genecode_data_fn | grep 'pseudogene' | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$4,$5,$9}' > $pseudogene_fn # this will actually report genes that have either gene_type or transcript_type that contains pseudogene word in it
gzip -f $pseudogene_fn
echo "Done getting genecode pseudogene into ${pseudogene_fn}.gz"

echo "Please open this code file if you want to see the commented out code that produces the count of different combinations of gene_type and transcript_type"
# zcat $genecode_data_fn | awk '{if (NR>5) print $0}' | awk -F'\t' -v cds="CDS" '{if ($3==cds) print $0}' | awk -F'\t' '{gsub("; ","\t",$9); print $0}' | awk -F'[\t ]' 'BEGIN { OFS="\t" } {gsub(/"/,"",$0); print $0}' | awk -F'[\t ]' 'BEGIN { OFS="\t" } { print $14,$18}' |sort|uniq -c

def check_gene_type_pseudo_gene(annot):
	annot = annot.split(';')
	if annot[-1] == '':
		annot = annot[:-1]
	annot = list(map(lambda x: tuple(x.strip().split()), annot))
	annot = dict(annot)
	return annot['gene_type']