import os 
import argparse
import pandas as pd
import numpy as np 
import glob
import pybedtools as bed
import helper
parser = argparse.ArgumentParser(description='Calculating the enrichment between the refernece chromatin state annotation and another segmentation states. This is useful when we want to calculate the overlap enrichment between full stack states and the states in other ct-spec models. Here, we will output the files in a way that replicates output of ChromHMM OverlapEnrichment')
parser.add_argument('--reference_segment_fn', type=str,
                    help='The segmentation fn that we will use to calculate the state overlap enrichment for. The rows of output correspond to states in this file')
parser.add_argument('--other_segment_fn', type=str,
                    help='other_segment_fn. The columns of output will be the states in this file')
parser.add_argument('--num_other_states', default=100, type=int,
                    help='number of chromHMM states in other_segment_fn')
parser.add_argument('--output_fn', type=str,
                    help = 'output_fn. Will look exactly like the output of ChromHMM OverlapEnrichment')
args = parser.parse_args()
print (args)
helper.check_file_exist(args.reference_segment_fn)
helper.check_file_exist(args.other_segment_fn)
helper.create_folder_for_file(args.output_fn)

def calculate_percentage_in_genome(reference_segment_fn):
	ref_segment_df = pd.read_csv(reference_segment_fn, header = None, index_col = None, sep = '\t')
	ref_segment_df.columns = ['chrom', 'start', 'end', 'state']
	result_df = pd.DataFrame()
	all_ref_states = np.unique(ref_segment_df['state'])
	all_ref_states = list(map(lambda x: int(x[1:]), all_ref_states)) # a list of unique state numbers
	all_ref_states = np.sort(all_ref_states) # 1--> 100, for full-stack states
	result_df['state (Emission order)'] = all_ref_states
	result_df.index = list(map(lambda x: 'E{}'.format(x), result_df['state (Emission order)']))
	ref_segment_df['length'] = ref_segment_df['end'] - ref_segment_df['start']
	gw_length = np.sum(ref_segment_df['length'])
	state_cover = ref_segment_df.groupby('state')['length'].sum()
	state_perc = state_cover/ gw_length * 100.0 # percentage of the genome that is in each state
	result_df['Genome %'] = state_perc
	result_df['num_bp_in_state'] = state_cover
	return result_df, gw_length

def calculate_overlap_enrichment_all_states(reference_segment_fn, other_segment_fn, num_other_states, output_fn):
	ref_segment_bed = bed.BedTool(reference_segment_fn)
	result_df, gw_length = calculate_percentage_in_genome(reference_segment_fn) # two columns so far: state (Emission order) and Genome %
	other_bed = bed.BedTool(other_segment_fn)
	inter_bed = ref_segment_bed.intersect(other_bed, wa = True, wb= True)
	inter_df = inter_bed.to_dataframe()
	inter_df.columns = ['ref_chrom', 'ref_start', 'ref_end', 'ref_state', 'other_chrom', 'other_start', 'other_end', 'other_state'] 
	inter_df['inter_length'] = inter_df.apply(lambda x: min(x['ref_end'], x['other_end']) - max(x['ref_start'], x['other_start']),axis = 1) # length of intersection between the reference and the other segmentation
	inter_df = inter_df.groupby(['ref_state', 'other_state'])['inter_length'].sum()
	inter_df = inter_df.to_frame().reset_index() # ref_state, other_state, inter_length
	inter_df = inter_df.pivot_table(values = 'inter_length', index = 'ref_state', columns = 'other_state') # columns: states in the other_bed, rows: states in ref_bed
	inter_df = inter_df.fillna(0)
	fract_gene_in_context = (inter_df.sum(axis = 0))/ gw_length # fraction of the genome that is in each of the contexts of the columns
	inter_df = inter_df.divide(result_df['num_bp_in_state'], axis = 0).divide(fract_gene_in_context, axis = 1) # FE = (#MS/#S) / (#M/#G)
	result_df = result_df.merge(inter_df, left_index = True, right_index = True) 
	result_df = result_df.drop('num_bp_in_state', axis = 1) # get rid of this column to make the output look like ChromHMM output
	result_df.loc[result_df.shape[0]] = ['Base', 100] + list(fract_gene_in_context*100.0) # the last row of result_df should show the percentage of the genome that is in each of the context
	result_df.to_csv(output_fn, header = True, index = False, sep = '\t')
	print('Done!')
	return 

calculate_overlap_enrichment_all_states(args.reference_segment_fn, args.other_segment_fn, args.num_other_states, args.output_fn)
