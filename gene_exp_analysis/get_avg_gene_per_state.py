import os 
import argparse
import pandas as pd
import numpy as np 
import glob
import pybedtools as bed
import helper
parser = argparse.ArgumentParser(description='Calculating avg gene expression in each state, in each available cell types. The gene_exp_folder should be where we store BingRen\'s gene expression data for mouse. segment_fn should be from ChromHMM model for the mouse')
parser.add_argument('--gene_exp_folder', type=str,
                    help='Where there gene exp data for different cell types are stored')
parser.add_argument('--segment_fn', type=str,
                    help='segment_fn. NOTE: Because our gene expression data is in mm9, we will need to pass the segment_fn in mm9 here')
parser.add_argument('--num_chromHMM_state', default=100, type=int,
                    help='number of chromHMM states')
parser.add_argument('--output_fn', type=str,
                    help = 'output_fn')
parser.add_argument('--state_annot_fn', type=str,
                    help = 'state_annot_fn')
args = parser.parse_args()
print (args)
helper.check_dir_exist(args.gene_exp_folder)
helper.check_file_exist(args.segment_fn)
helper.create_folder_for_file(args.output_fn)
helper.check_file_exist(args.state_annot_fn)
SEGMENT_LENGTH = 200

def open_gene_exp_data_one_ct(gene_exp_fn):
	df = pd.read_csv(gene_exp_fn, header = 0, index_col = None, sep = '\t')
	df = df[['chr', 'left', 'right', 'FPKM']]
	df['FPKM'] = df['FPKM'].apply(lambda x: np.log(x+1))
	df.columns = ['chrom', 'start', 'end', 'logFPKM']
	return df

def get_avg_gene_exp_one_ct(segment_fn, gene_exp_fn, num_chromHMM_state):
	segment_bed = bed.BedTool(segment_fn)
	exp_df = open_gene_exp_data_one_ct(gene_exp_fn)
	exp_bed = bed.BedTool.from_dataframe(exp_df)
	segment_bed = segment_bed.intersect(exp_bed, wa = True, wb = True)
	comb_df = segment_bed.to_dataframe()
	comb_df.columns = ['sChrom', 'sStart', 'sEnd', 'state', 'eChrom', 'eStart', 'eEnd', 'logFPKM']
	comb_df['gene_length'] = comb_df['eEnd'] - comb_df['eStart']
	comb_df['gene_length_segments'] = comb_df['gene_length'] / SEGMENT_LENGTH
	comb_df['com_start'] = comb_df[['sStart', 'eStart']].max(axis = 1) #start of the region we are interested is the start of either the gene or the segment of chromatin state. Whichever is greater is the start of the intersection. 
	comb_df['com_end'] = comb_df[['sEnd', 'eEnd']].min(axis = 1) # same argument as above with the end of the region 
	comb_df['num_segments'] = (comb_df['com_end'] - comb_df['com_start']) / SEGMENT_LENGTH # number of segments that are shared between the genes and the chromatin state
	comb_df['num_segments'] = comb_df['num_segments'].apply(np.ceil) # round up the number of segments
	# Number of segments is the number of basepair in the region we are interested in / number of basepairs per segment bin
	comb_df = comb_df[[u'state', u'gene_length', u'com_start', u'com_end', u'num_segments', 'gene_length_segments', 'logFPKM']] 
	# now onto processing the gene expression data
	result_avg_exp_S = pd.Series([], dtype = float)
	for state_index in range(num_chromHMM_state):
		one_based_state_index = state_index + 1
		this_state_org_df = (comb_df[comb_df['state'] == 'E' + str(one_based_state_index)])
		bp_unif_exp_S = this_state_org_df['logFPKM'] * this_state_org_df['num_segments'] / this_state_org_df['gene_length_segments'] # a pandas series where each entry correspond to a position on the genome annotated as this state
		total_weighted_segments = np.sum(this_state_org_df['num_segments'] / this_state_org_df['gene_length_segments']) # a number
		avg_exp_bp_unif = bp_unif_exp_S.sum() / total_weighted_segments
		result_avg_exp_S['E' + str(one_based_state_index)] = avg_exp_bp_unif
	return result_avg_exp_S # a pandas series with index: E1--> E100, values: avg gene expression in each state in the current cell type

def get_cell_name_from_fn(fn):
	fn = fn.split('/')[-1].split('.expr')[0].split('.gene')[0]
	first_fn = fn.split('-')[0]
	if first_fn.endswith('1') or first_fn.endswith('2') :
		return first_fn
	last_fn = fn.split('-')[-1]
	if (last_fn == '1' or last_fn == '2' or last_fn == '3') and first_fn != 'heart' and first_fn != 'liver':
		return first_fn + last_fn
	if (first_fn == 'heart' or first_fn == 'liver'):
		return '_'.join(fn.split('-'))
	else:
		return '_'.join(fn.split('-')[:2])
	return ''

def get_avg_gene_exp_all_ct(segment_fn, gene_exp_folder, num_chromHMM_state, output_fn):
	gene_exp_fn_list = glob.glob(gene_exp_folder+'/*.expr')
	ct_list = list(map(get_cell_name_from_fn, gene_exp_fn_list))
	result_index = list(map(lambda x: 'E{}'.format(x+1), range(num_chromHMM_state)))
	result_avg_exp_df = pd.DataFrame(columns = ct_list, index = result_index)
	for ct_index, ct in enumerate(ct_list):
		result_avg_exp_df[ct] = get_avg_gene_exp_one_ct(segment_fn, gene_exp_fn_list[ct_index], num_chromHMM_state)
	result_avg_exp_df.to_csv(output_fn, header = True, index = True, sep = '\t', compression = 'gzip')
	print('Done!')
	return 

get_avg_gene_exp_all_ct(args.segment_fn, args.gene_exp_folder, args.num_chromHMM_state, args.output_fn)