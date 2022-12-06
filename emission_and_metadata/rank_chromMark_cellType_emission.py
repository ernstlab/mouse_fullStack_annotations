import pandas as pd
import numpy as np
import sys
import os
import helper

def prepare_color_dict():
	assay_meta_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/emissions/assay_count.csv'
	organ_meta_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/emissions/organ.txt'
	assay_df = pd.read_csv(assay_meta_fn, header = 0, index_col = None, sep = ',')[['mark', 'color', 'big_group']]
	assay_df = assay_df.rename(columns = {'color': 'mark_color', 'big_group': 'assay_big_group'})
	assay_color_dict = dict(assay_df.groupby(['mark', 'mark_color']).groups.keys()) # keys: assay, values: colors
	organ_df = pd.read_csv(organ_meta_fn, header = None, index_col = None, sep = '\t')
	organ_df.columns = ['group', 'color', 'order']
	organ_df['group'] = organ_df['group'].apply(lambda x: x.lower())
	organ_color_dict = dict(organ_df.groupby(['group', 'color']).groups.keys()) # keys: biomsaple_name, values: color
	return assay_color_dict, organ_color_dict

ASSAY_COLOR_DICT, organ_COLOR_DICT = prepare_color_dict()

def get_ranked_mark_name (row_data):
	sorted_rank = row_data.sort_values(ascending = False) # 1 --> n
	return pd.Series(sorted_rank.index) # get the index whcih is the list of experiments, ordered from least emitted to most emitted in each state. Convert to pd.Series so that we can concatenate into data frame later.

def get_ranked_exp_df (fn) : # fn should be emission_fn
	df = pd.read_csv(fn, header = 0, sep = '\t') # get emission df
	df = df.rename(columns = {'State (Emission order)' : 'state'}) # one column is state, others  are all experiment names
	df.index = df['state'] # state column becomes the index of the dataframe
	df = df[df.columns[1:]] # get rid of state column so that we can rank the emission probabilities of each state
	rank_df = df.rank(axis = 1) # index: state, columns: experiments ordered just like in df. Cells: rank of each experiment within each state, 1: least emission --> n: highest emission
	rank_df = rank_df.apply(get_ranked_mark_name, axis = 1) # index; state, columns: rank 1 --> n most emitted experiment to least emitted experiments
	rank_df.columns = map(lambda x: 'r_' + str(x + 1), range(len(rank_df.columns)))
	return rank_df


def color_organ_names(val):
	if val == "":
		color = organ_COLOR_DICT["NA"]
	else:
		color = organ_COLOR_DICT[val]
	return 'background-color: %s' % color

def color_mark_names(val):
	if val == "":
		color = ASSAY_COLOR_DICT['NaN']
	else:
		color = ASSAY_COLOR_DICT[val]
	return 'background-color: %s' % color

def read_exp_organ_dict():
	meta_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/emissions/metadata_from_jason.txt'
	meta_df = pd.read_csv(meta_fn, header = 0, index_col = None, sep = '\t') # 'experimentID', 'organ_slims', 'cell_slims', 'developmental_slims', 'system_slims', 'biosample_summary', 'simple_biosample_summary'
	meta_df['organ_group'] = meta_df['organ_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('musculature_of_body', 'musculature').replace('', 'unknown')
	return dict(zip(meta_df.experimentID, meta_df.organ_group))

def get_painted_excel_ranked_exp(rank_df, output_fn, num_top_marks, organ_color_dict, assay_color_dict):
	EXP_ORGAN_DICT = read_exp_organ_dict()
	organ_df = rank_df.applymap(lambda x: EXP_ORGAN_DICT[x.split('_')[2]]) # first convert the data to only contain the experimentID, then convert from experimentID to ORGAN
	organ_df.reset_index(inplace = True)
	chrom_mark_df = rank_df.applymap(lambda x: x.split('_')[-2].split('-')[0]) # H3K9me3
	chrom_mark_df.reset_index(inplace = True)
	columns_to_paint = list(map(lambda x: 'r_' + str(x + 1), range(num_top_marks)))
	organ_df = organ_df[['state'] + columns_to_paint]
	chrom_mark_df = chrom_mark_df[['state'] + columns_to_paint]
	colored_organ_df = organ_df.style.applymap(color_organ_names, subset = columns_to_paint)
	colored_chrom_mark_df = chrom_mark_df.style.applymap(color_mark_names, subset = columns_to_paint)
	# save file
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	colored_organ_df.to_excel(writer, sheet_name = 'cell_group')
	colored_chrom_mark_df.to_excel(writer, sheet_name = 'chrom_mark')
	writer.save()
	print ("Done saving data into " + output_fn)



def main():
	if len(sys.argv) != 4:
		usage()
	emission_fn = sys.argv[1]
	output_fn = sys.argv[2]
	helper.create_folder_for_file(output_fn)
	num_top_marks = helper.get_command_line_integer(sys.argv[3])
	print ('Done getting command line argument')
	assay_color_dict, organ_color_dict = prepare_color_dict()
	rank_df = get_ranked_exp_df(emission_fn) # index: states, columns: rank of experiments r_1: highest emitted, r_n: lowest emitted --> cells: names of experiments. Example: E118-H3K9me3
	print ("Done getting ranked data")
	get_painted_excel_ranked_exp(rank_df, output_fn, num_top_marks, organ_color_dict, assay_color_dict)



def usage():
	print ("python rank_chrom_mark_celltype_emission.py")
	print ("emission_fn")
	print ("output_fn: execl file where the output data of ranked cell type and chrom marks are stored")
	print ("num_top_marks: number of top marks that we want to report")
	exit(1)

main()
