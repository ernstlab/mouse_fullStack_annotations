import os
import pandas as pd 
import helper
import numpy as np 
import pybedtools as bed 
import argparse
import seaborn as sns

parser = argparse.ArgumentParser(description="Find positive pairs of variants")
parser.add_argument("--segment_fn", type=str, required=True, help="full-stack segmentation file")
parser.add_argument('--lecif_fn', type = str, required = True, help='bed file of lecif score')
parser.add_argument("--state_annot_fn", type=str, required=True, help="File of state characterization")
parser.add_argument("--output_folder", type = str, required=True, help = 'output_folder')
parser.add_argument('--num_state', type = int, required=False, default=100, help='number of states in the segmentation file')

def get_state_annot_df(state_annot_fn):
    # state_annot_fn = '../../..//ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv'
    try:
        state_annot_df = pd.read_csv(state_annot_fn, sep = ',', header = 0)
    except:
        state_annot_df = pd.read_csv(state_annot_fn, sep = '\t', header = 0)
    state_annot_df = state_annot_df[['state', 'color', 'mneumonics', 'state_order_by_group', 'comments']]
    return state_annot_df

def color_state_annotation(row_data, index_to_color_list):
    results = [""] * len(row_data) # current paint all the cells in the rows with nothing, no format yet ---> all white for now
    state_annot_color = row_data['color']
    if not pd.isna(row_data['color']):
        for index in index_to_color_list:
            results[index] = 'background-color: %s' % state_annot_color # the third cell from the left is the state annotation cells
    return results

def transform_state_to_number(state):
	'''
	E1 --> 1
	'''
	if state.startswith('E'):
		return int(state[1:])
	return int(state)
	

def combine_lecif_df_with_state_annot_df(lecif_df, state_annot_fn):
    state_annot_df = get_state_annot_df(state_annot_fn)
    lecif_df['state'] = (lecif_df['state']).apply(transform_state_to_number)
    lecif_df = lecif_df.merge(state_annot_df, how = 'left', left_on = 'state', right_on = 'state', suffixes = ('_x', '_y'))
    lecif_df = lecif_df.sort_values(by = 'state_order_by_group')
    lecif_df = lecif_df.reset_index()  
    lecif_df = lecif_df.drop(columns = ['index'])
    lecif_df['mneumonics'] = lecif_df['mneumonics'].astype(str)
    lecif_df['color'] = lecif_df['color'].astype(str)
    lecif_df['state_order_by_group'] = lecif_df['state_order_by_group'].astype(int)
    return lecif_df

def get_enrichment_colored_df(lecif_df, save_fn, state_annot_fn):
    cm = sns.light_palette("blue", as_cmap=True)
    lecif_df = combine_lecif_df_with_state_annot_df(lecif_df, state_annot_fn)
    mneumonics_index_in_row = lecif_df.columns.get_loc('mneumonics') # column index, zero-based of the menumonics entries, which we will use to paint the right column later
    colored_df = lecif_df.style.background_gradient(subset = pd.IndexSlice[:, ['mean']], cmap = cm)
    colored_df = colored_df.apply(lambda x: color_state_annotation(x, [mneumonics_index_in_row]), axis = 1)
    colored_df.to_excel(save_fn, engine = 'openpyxl')
    return colored_df


def calculate_avg_lecif_one_state(segment_df, lecif_bed, state):
	'''
	state should be E1 --> E100
	'''
	state_bed = segment_df[segment_df[3] == state]
	state_bed = bed.BedTool.from_dataframe(state_bed)
	intersect_bed = state_bed.intersect(lecif_bed, wa = True, wb = True)
	intersect_bed = intersect_bed.to_dataframe()
	intersect_bed.columns = ['chrom_s', 'start_s', 'end_s', 'state', 'chrom_lecif', 'start_lecif', 'end_lecif', 'lecif']
	result = intersect_bed['lecif'].describe()
	result['state'] = state
	return result # mean, std. min, 25%, 50%, 75%, max

def get_avg_lecif(segment_fn, lecif_fn, num_state, output_folder):
	segment_df = pd.read_csv(segment_fn, header = None, index_col = None, sep = '\t')
	lecif_bed = bed.BedTool(lecif_fn)
	result_df = pd.DataFrame(columns =['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max', 'state'])
	for state in range(num_state):
		state = 'E{}'.format(state+1)
		state_stats = calculate_avg_lecif_one_state(segment_df, lecif_bed, state)
		result_df.loc[result_df.shape[0]] = state_stats
	result_df.drop('count', axis = 1, inplace = True)
	output_fn = os.path.join(output_folder, 'avg_lecif_scores.txt')
	result_df.to_csv(output_fn, header = True, index = False, sep = '\t')
	return result_df

if __name__ == '__main__':
	args = parser.parse_args()
	helper.check_file_exist(args.segment_fn)
	helper.check_file_exist(args.lecif_fn)
	helper.make_dir(args.output_folder)
	helper.check_file_exist(args.state_annot_fn)
	avg_lecif_fn = os.path.join(args.output_folder, 'avg_lecif_scores.txt')
	if not os.path.isfile(avg_lecif_fn):
		avg_lecif_df = get_avg_lecif(args.segment_fn, args.lecif_fn, args.num_state, args.output_folder)
	else:
		avg_lecif_df = pd.read_csv(avg_lecif_fn, header = 0, index_col = None, sep = '\t')
	excel_save_fn = os.path.join(args.output_folder, 'avg_lecif_scores.xlsx')
	get_enrichment_colored_df(avg_lecif_df, excel_save_fn, args.state_annot_fn)	
	print ('Done!')

