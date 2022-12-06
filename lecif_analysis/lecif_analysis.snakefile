import numpy as np 
lecif_raw_fn = '/u/home/h/havu73/project-ernst/data/LECIF/mm10/mm10.LECIFv1.1.bed.gz'
bed_folder = '/u/home/h/havu73/project-ernst/data/LECIF/mm10/bed_range/'
bin_width = 0.1
upper_bound_list = np.around(np.arange(bin_width, 1.1, bin_width), decimals = 1)
rule all:
	input:
		expand(os.path.join(bed_folder, 'lecif_{ub}.bed.gz'), ub = upper_bound_list),

rule get_lecif_bed_one_range:
	input:
		lecif_raw_fn,
	output:
		os.path.join(bed_folder, 'lecif_{ub}.bed.gz'),
	params:
		output_fn_no_gz = os.path.join(bed_folder, 'lecif_{ub}.bed')
	shell:
		"""
		lb=$( echo "{wildcards.ub}-{bin_width}" |bc )
		echo $lb
		zcat {lecif_raw_fn} | awk -F'\t' -v lower_bound=$lb '{{if ($4>=lower_bound && $4<{wildcards.ub}) print $0}}' > {params.output_fn_no_gz}
		gzip {params.output_fn_no_gz}
		"""