import os
import glob
full_stack_segment_fn = 'genome_100_segments.bed.gz' # should be replaced with where you store the 100-state chromatin state annotation
ct_segment_folder = './example_perCT_segment/' # should be replaced with the folder path to where you store the per_ct annotations. NOTE: the code assume that your per-ct annotation (bed files) are sorted, if not you need to use the command line: zcat <bed_file_name> | sort -k1,1 k2,2n > <output_fn> 
ct_list = glob.glob(ct_segment_folder + '/*_15_segments.bed.gz') # we assume that within the folder ct_segment_folder, the files are named based on format <ct_name>_15_segments.bed.gz
ct_list = list(map(lambda x: x.split('/')[-1].split('_15_segments.bed.gz')[0], ct_list)) 
NUM_SAMPLE_SEGMENT_PER_STATE = 100
output_folder = './output//random_represent_with_15state' # output folder of analysis of the estimated overlap probability for full-stack states with per-ct states
overlap_output_folder = './output//overlap_with_ct_segment' # output folder of the analysis of overlap enrichment of full-stack states with per-ct states
ct_state_annot_fn = './example_perCT_segment/15_state_annotations.csv' # where we can get the metadata of the states within the per-ct chromatin state model
full_state_annot_fn = '../state_annotation_processed.csv' # path to the file of full-stack states' metadata (state characterizations)
ct_group_fn = './example_perCT_segment/tissue_annotation.csv'
num_state_ct_model = 15

seed_list = [800, 922,  23, 204, 132, 992,  60, 650, 761, 154, 432, 760, 999, 969, 955, 986, 981, 246,  35, 438, 116]
rule all:
	input:
		os.path.join(output_folder, 'summary', 'avg_prop_stateCT_per_group.xlsx'),
		os.path.join(overlap_output_folder, 'summary', 'max_enriched_states.xlsx'),

rule sample_segment_equal_per_state:
	# sample regions on the genome such that for each state, the number of sampled regions are similar
	input:
		full_stack_segment_fn, # from raw data
	output:
		(os.path.join(output_folder, 'seed_{seed}', 'temp_sample_region.bed.gz'))
	shell:
		"""
		python sample_region_for_state_representation.py {input} {NUM_SAMPLE_SEGMENT_PER_STATE} {output} {wildcards.seed}
		"""

rule overlap_sample_region_with_ct_segment:
	input:	
		(os.path.join(output_folder, 'seed_{seed}', 'temp_sample_region.bed.gz')), # from rule sample_segment_equal_per_state
		expand(os.path.join(ct_segment_folder, '{ct}_15_segments.bed.gz'), ct = ct_list) # Note we made the assumption that these files are all sorted 'sort -k1,1 -k2,2n'
	output:
		os.path.join(output_folder, 'seed_{seed}', 'sample_segment_fullStack_ctState.bed.gz'),
	params:
		ct_list_string = " ".join(ct_list),
		output_no_gz = os.path.join(output_folder, 'seed_{seed}', 'sample_segment_fullStack_ctState.bed')
	shell:
		"""
		command="zcat {input[0]} | sort -k1,1 -k2,2n " # first sort the file of sample regions, where each state appear NUM_SAMPLE_SEGMENT_PER_STATE times
		output_header="chrom\\tstart\\tend\\tfull_stack"
		rm -f {output} # so we can overwrite it
		rm -f {params.output_no_gz} # so we can overwrite it
		for ct in {params.ct_list_string}
		do 
			ct_segment_fn={ct_segment_folder}/${{ct}}_15_segments.bed.gz
			command="$command | bedtools map -a stdin -b ${{ct_segment_fn}} -c 4 -o collapse"
			output_header="${{output_header}}\\t${{ct}}"
		done
		command="$command >> {params.output_no_gz} "
		# now onto the writing part
		echo  -e $output_header > {params.output_no_gz} # write the header first # need the -e flag for tabs to be considered seriously
		eval $command # write the content
		gzip {params.output_no_gz} # gzip
		""" 

rule calculate_overlap_enrichment_with_ct_segment:
	input:
		full_stack_segment_fn,
		os.path.join(ct_segment_folder, '{ct}_15_segments.bed.gz'),
	output:
		os.path.join(overlap_output_folder, '{ct}', 'overlap_enrichment.txt'),
	params:
		code='/Users/vuthaiha/Desktop/window_hoff/source/25state_enrichments/overlap_enrichment_with_ct_spec_state.py',
		num_other_states = 15,
	shell:
		"""
		python {params.code} --reference_segment_fn {full_stack_segment_fn} --other_segment_fn {input[1]} --num_other_states {params.num_other_states} --output_fn {output}
		"""

rule calculate_summary_overlap_with_ct_segment:
	input:
		expand(os.path.join(overlap_output_folder, '{ct}', 'overlap_enrichment.txt'), ct = ct_list),
	output:
		os.path.join(overlap_output_folder, 'summary', 'max_enriched_states.xlsx'),
	shell:
		"""
		python /Users/vuthaiha/Desktop/window_hoff/source/mm10_annotations/relationship_with_ct_spec_state/get_ranked_max_enriched_25_state.py --input_folder {overlap_output_folder} --output_fn {output} --num_state_ct_model {num_state_ct_model} --full_state_annot_fn {full_state_annot_fn} --ct_state_annot_fn {ct_state_annot_fn} --ct_group_fn {ct_group_fn}
		"""

rule calculate_summary_prob_overlap_with_ct_segment:
	input:
		expand(os.path.join(output_folder, 'seed_{seed}', 'sample_segment_fullStack_ctState.bed.gz'), seed = seed_list)
	output:
		os.path.join(output_folder, 'summary', 'avg_prop_stateCT_per_group.xlsx')
	params:
		prob_summ_folder = os.path.join(output_folder, 'summary')
	shell:
		"""
		python /Users/vuthaiha/Desktop/window_hoff/source/mm10_annotations/relationship_with_ct_spec_state/calculate_summary_sample_regions.py --all_seed_folder {output_folder} --output_folder {params.prob_summ_folder} --num_state_ct_model {num_state_ct_model} --full_state_annot_fn {full_state_annot_fn} --ct_state_annot_fn {ct_state_annot_fn} --ct_group_fn {ct_group_fn}
		"""