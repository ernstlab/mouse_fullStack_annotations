import pandas as pd
import seaborn as sns
import numpy as np
import os
import sys
import helper
import glob
import argparse
parser = argparse.ArgumentParser(description = 'Create an excel of the top ct-spec states enriched in each full-stack state')
parser.add_argument('--input_folder', type=str,
	help = "where there are multiple subfolders, each containing enrichment data for different cell type specific model")
parser.add_argument('--output_fn', type=str,
	help = 'Where the ct-state with maximum fold enrichment across all ct-state states are reported for all cell types for each state is reported')
parser.add_argument('--num_state_ct_model', type=int,
	help = 'number of states in the ct-spec models')
parser.add_argument('--full_state_annot_fn', type=str,
	help = 'file of the full-stack state annotations')
parser.add_argument('--ct_state_annot_fn', type=str,
	help = 'file of the ct-spec state annotations')
parser.add_argument('--ct_group_fn', type=str,
	help = 'file of the annotations of the cell types')
args = parser.parse_args()
print(args)
helper.check_dir_exist(args.input_folder)
helper.create_folder_for_file(args.output_fn)
helper.check_file_exist(args.full_state_annot_fn)
helper.check_file_exist(args.ct_state_annot_fn)
helper.check_file_exist(args.ct_group_fn)

def get_full_state_annot_df(full_state_annot_fn):
	df = pd.read_csv(full_state_annot_fn, header = 0, index_col = None, sep = '\t')
	FULL_STATE_COLOR_DICT = pd.Series(df.color.values, index = df.mneumonics)
	df = df[['state', 'mneumonics', 'state_order_by_group']]
	return FULL_STATE_COLOR_DICT, df

FULL_STATE_COLOR_DICT, FULL_STATE_ANNOT_DF = get_full_state_annot_df(args.full_state_annot_fn) #FULL_STATE_ANNOT_DF: state, mneumonics, state_order_by_group, FULL_STATE_COLOR_DICT: keys mneumonics, values color

def get_celltype_color_map(ct_group_fn):
	df = pd.read_csv(ct_group_fn, header = 0, index_col = None, sep = '\t')
	CELLTYPE_COLOR_MAP = pd.Series(df.color.values, index = df.tissue_stage).to_dict() # convert two columns in to a dictionary: keys: cell type, values: color corresponding to that tissue
	CELLTYPE_COLOR_MAP['NA'] = '#E5B8E8'
	CELLTYPE_CELL_GROUP_MAP = pd.Series(df.group.values, index = df.tissue_stage).to_dict() # keys: cell type name, values: cell group
	group_color_df = (df[['group', 'color']]).drop_duplicates()
	CELL_GROUP_COLOR_CODE_DICT = pd.Series(group_color_df.color.values, index = group_color_df.group).to_dict() # keys: group, values: color
	return CELLTYPE_COLOR_MAP, CELLTYPE_CELL_GROUP_MAP, CELL_GROUP_COLOR_CODE_DICT, df

CELLTYPE_COLOR_MAP, CELLTYPE_CELL_GROUP_MAP, CELL_GROUP_COLOR_CODE_DICT, CT_DF = get_celltype_color_map(args.ct_group_fn)

def get_ct_spec_state_annot(ct_state_annot_fn):
	df = pd.read_csv(ct_state_annot_fn, header = 0, index_col = None, sep = '\t') # state, mneumonics, meaning, color
	CT_STATE_MNEUNOMICS_LIST = list(df.mneumonics)
	CT_STATE_COLOR_DICT = pd.Series(df.color.values, index = df.mneumonics).to_dict() # keys: state mneumonics, values: state color
	return CT_STATE_MNEUNOMICS_LIST, CT_STATE_COLOR_DICT

CT_STATE_MNEUNOMICS_LIST, CT_STATE_COLOR_DICT = get_ct_spec_state_annot(args.ct_state_annot_fn)

def get_rid_of_stupid_file_tail(context_name):
	if context_name.endswith('.bed.gz'):
		return(context_name[:-7])
	else:
		return(context_name)

def color_full_state_names(val):
	if val == "":
		color = '#FFFFFF'
	else:
		color = FULL_STATE_COLOR_DICT[val]
	return 'background-color: %s' % color

def color_cell_group_names(val):
	if val == "":
		color = '#FFFFFF'
	else:
		color = CELLTYPE_COLOR_MAP[val]
	return 'background-color: %s' % color

def color_tissue_names(val):
	if val == "":
		color = '#FFFFFF'
	else:
		color = CELL_GROUP_COLOR_CODE_DICT[val]
	return 'background-color: %s' % color

def get_one_enrichment_ct_model_df(fn, num_state_ct_model):
	df = pd.read_csv(fn, sep = '\t', header = 0)
	df = df.rename(columns = {'state (Emission order)' : 'state', 'Genome %': 'percent_in_genome'}) # rename some columns so that it is easier to write column names later
	state_colName_list = list(map(lambda x: 'U' + str(x + 1), range(num_state_ct_model))) #CUSTOM: this is a line of code that could be improved for better customization
	df = df[['state', 'percent_in_genome'] + state_colName_list] # get the data frame to display columns in the expected order
	df.columns = ['state', 'percent_in_genome'] + CT_STATE_MNEUNOMICS_LIST # rename the columns so that instead of 'state1.bed.gz' we have '1_TssA'
	df['max_enrichment'] = (df.drop(['state', 'percent_in_genome'], axis = 1)).max(axis = 1) # find the maximum enrichment values in this row 
	df['max_enrichment_context'] = (df.drop(['state', 'percent_in_genome'], axis = 1)).idxmax(axis = 1)
	(nrow, ncol) = df.shape
	df = df.drop(nrow - 1) # drop the last row, which is the 'Base' row with percentage that each enrichment context occupies the genome
	return df

def get_25_state_annot(state):
	# 'state_1' --> ''1_TssA''
	state_index = int(state.split('_')[1]) - 1 # zero-based
	return CT_STATE_MNEUNOMICS_LIST[state_index]

def get_CT_STATE_COLOR_DICT(state):
	# '1_TssA' --> 'red'
	color = CT_STATE_COLOR_DICT[state]
	return 'background-color: %s' % color

def get_max_enrichment_25_state_df_ordered_time_stamp(max_enrich_25_state_all_ct_df):
	"""
	max_enrich_25_state_all_ct_df: rows: full stack states, columns: state of the 25-state system that is most enriched with the full-stack states WITHIN a cell type. Cell type are the column names  
	--> a plot where cell types are juxtaposed based on the cell groups they belong to. The output dataframe will be colored properly
	"""
	cell_group_df = pd.DataFrame(columns = ['tissue_stage']) # rows: the cell types that are used in this analysis (colnames of max_enrich_25_state_all_ct_df), columns : 
	cell_group_df['tissue_stage'] = max_enrich_25_state_all_ct_df.columns[1:] # we skip the first column because that's 'state'
	cell_group_df = pd.merge(cell_group_df, CT_DF, how = 'left', left_on = 'tissue_stage', right_on = 'tissue_stage') # merge the cell types so that we can get the information that we want
	# Now let's get the count of the number of cell types that are of each cell groups, and then get the unique cell groups, ordered by descending counts
	unique_cell_groups = ['forebrain', 'midbrain', 'hindbrain', 'neural-tube', 'facial-prominence', 'limb', 'intestine', 'stomach', 'liver', 'kidney', 'lung', 'heart']#CUSTOM: change this line of code for better customization	
	unique_stages = ['P0', 'e10.5', 'e11.5', 'e12.5', 'e13.5', 'e14.5', 'e15.5', 'e16.5']#CUSTOM: change this line of code for better customization
	cell_group_df['ct'] = pd.Categorical(cell_group_df['ct'], unique_cell_groups) # so that we can easily sort them later based on this custom order
	# list of cell types, arranged such that those of the same cell_groups are juxtaposed
	rearranged_tissue_stages = []
	for stage in unique_stages:
		this_cg_df = cell_group_df[cell_group_df['stage'] == stage] # filter out rows with this cellgroup
		this_cg_df = this_cg_df.sort_values(['ct'], ignore_index = True)
		rearranged_tissue_stages += list(this_cg_df['tissue_stage'])
	# now we will rearrange the columns of max_enrich_25_state_all_ct_df based on the order that we just got from rearranged_tissue_stages
	max_enrich_25_state_all_ct_df['state'] = max_enrich_25_state_all_ct_df['state'].astype(int)
	max_enrich_25_state_all_ct_df = max_enrich_25_state_all_ct_df.merge(FULL_STATE_ANNOT_DF, how = 'left', left_on = 'state', right_on = 'state')
	max_enrich_25_state_all_ct_df = max_enrich_25_state_all_ct_df[['state'] + rearranged_tissue_stages + ['mneumonics', 'state_order_by_group']]
	max_enrich_25_state_all_ct_df = max_enrich_25_state_all_ct_df.sort_values('state_order_by_group', ignore_index = True)
	last_row_index = max_enrich_25_state_all_ct_df.shape[0]
	max_enrich_25_state_all_ct_df.loc[last_row_index] = [''] + rearranged_tissue_stages + ['', ''] # add one more row that specific the cell_group of each of the cell types
	colored_25_state_df = max_enrich_25_state_all_ct_df.style.applymap(color_cell_group_names, subset = pd.IndexSlice[last_row_index, :]) # color the row that contain the cell_group of cell types
	colored_25_state_df = colored_25_state_df.applymap(get_CT_STATE_COLOR_DICT, subset = pd.IndexSlice[:(last_row_index - 1) , colored_25_state_df.columns[1:max_enrich_25_state_all_ct_df.shape[1]-2]]) # color the states of the 25-system that is most enriched in each of the full-stack states, exclude the first column and the last two columns
	colored_25_state_df = colored_25_state_df.applymap(color_full_state_names, subset = pd.IndexSlice[:,['mneumonics']]) # color full stack state names
	return colored_25_state_df


def get_max_enrichment_25_state_df_ordered_cell_group(max_enrich_25_state_all_ct_df):
	"""
	max_enrich_25_state_all_ct_df: rows: full stack states, columns: state of the 25-state system that is most enriched with the full-stack states WITHIN a cell type. Cell type are the column names  
	--> a plot where cell types are juxtaposed based on the cell groups they belong to. The output dataframe will be colored properly
	"""
	cell_group_df = pd.DataFrame(columns = ['tissue_stage']) # rows: the cell types that are used in this analysis (colnames of max_enrich_25_state_all_ct_df), columns : 
	cell_group_df['tissue_stage'] = max_enrich_25_state_all_ct_df.columns[1:] # we skip the first column because that's 'state'
	cell_group_df = pd.merge(cell_group_df, CT_DF, how = 'left', left_on = 'tissue_stage', right_on = 'tissue_stage') # merge the cell types so that we can get the information that we want
	# Now let's get the count of the number of cell types that are of each cell groups, and then get the unique cell groups, ordered by descending counts
	unique_cell_groups = ['forebrain', 'midbrain', 'hindbrain', 'neural-tube', 'facial-prominence', 'limb', 'intestine', 'stomach', 'liver', 'kidney', 'lung', 'heart']#CUSTOM: change this line of code for better customization
	unique_stages = [ 'P0', 'e10.5', 'e11.5', 'e12.5', 'e13.5', 'e14.5', 'e15.5', 'e16.5']#CUSTOM: change this line of code for better customization
	cell_group_df['stage'] = pd.Categorical(cell_group_df['stage'], unique_stages) # so that we can easily sort them later based on this custom order
	# list of cell types, arranged such that those of the same cell_groups are juxtaposed
	rearranged_tissue_stages = []
	for cg in unique_cell_groups:
		this_cg_df = cell_group_df[cell_group_df['ct'] == cg] # filter out rows with this cellgroup
		this_cg_df = this_cg_df.sort_values(['stage'], ignore_index = True)
		rearranged_tissue_stages += list(this_cg_df['tissue_stage'])
	max_enrich_25_state_all_ct_df['state'] = max_enrich_25_state_all_ct_df['state'].astype(int) # so that we can merge with the df of full stack state annotation
	# now we will rearrange the columns of max_enrich_25_state_all_ct_df based on the order that we just got from rearranged_tissue_stage
	max_enrich_25_state_all_ct_df = max_enrich_25_state_all_ct_df.merge(FULL_STATE_ANNOT_DF, how = 'left', left_on = 'state', right_on = 'state')
	max_enrich_25_state_all_ct_df = max_enrich_25_state_all_ct_df[['state'] + rearranged_tissue_stages + ['mneumonics', 'state_order_by_group']]
	max_enrich_25_state_all_ct_df = max_enrich_25_state_all_ct_df.sort_values('state_order_by_group', ignore_index = True)
	last_row_index = max_enrich_25_state_all_ct_df.shape[0]
	max_enrich_25_state_all_ct_df.loc[last_row_index] = [''] + rearranged_tissue_stages + ['', ''] # add one more row that specific the cell_group of each of the cell types
	colored_25_state_df = max_enrich_25_state_all_ct_df.style.applymap(color_cell_group_names, subset = pd.IndexSlice[last_row_index, :]) # color the row that contain the cell_group of cell types
	colored_25_state_df = colored_25_state_df.applymap(get_CT_STATE_COLOR_DICT, subset = pd.IndexSlice[:(last_row_index - 1) , colored_25_state_df.columns[1:max_enrich_25_state_all_ct_df.shape[1]-2]]) # color the states of the 25-system that is most enriched in each of the full-stack states, exclude the first column and the last two columns
	colored_25_state_df = colored_25_state_df.applymap(color_full_state_names, subset = pd.IndexSlice[:,['mneumonics']]) # color full stack state names
	return colored_25_state_df


def get_rank_25_state_df(max_enrich_25_state_all_ct_df, max_enrich_value_all_ct_df, output_fn):
	(num_full_state, num_ct) = max_enrich_25_state_all_ct_df.shape # nrow, ncol
	# color 25_state_cell_group_df: df where cell types of the same cell group are juxtaposed
	colored_25_state_cell_group_df = get_max_enrichment_25_state_df_ordered_cell_group(max_enrich_25_state_all_ct_df)
	# color colored_25_state_stage_df: df where the samples of the same time stamps are juxtaposed 
	colored_25_state_stage_df = get_max_enrichment_25_state_df_ordered_time_stamp(max_enrich_25_state_all_ct_df)
	# save the excel
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	colored_25_state_cell_group_df.to_excel(writer, sheet_name = 'cell_group_25_state')
	colored_25_state_stage_df.to_excel(writer, sheet_name = 'stage_25_state')
	writer.save()



def get_ranked_max_enriched_25_state(input_folder, output_fn, num_state_ct_model):
	all_ct_fn_list = glob.glob(input_folder + '/*/overlap_enrichment.txt')
	all_ct_list = list(map(lambda x: x.split('/')[-2], all_ct_fn_list)) # from '/path/to/E129/overlap_enrichment.txt' to 'E129' 
	all_ct_df_list = list(map(lambda x: get_one_enrichment_ct_model_df(x, num_state_ct_model), all_ct_fn_list))

	max_enrich_value_all_ct_df = pd.DataFrame({'state' : (all_ct_df_list[0])['state']}) # create a data frame with only a column called state that are the states in the full-stack model. This data frame store the enrichment values
	max_enrich_25_state_all_ct_df = pd.DataFrame({'state' : (all_ct_df_list[0])['state']}) # this data frame stores the names of states that are most enriched with each of the full-stack state in each of the cell type that we look at
	for ct_index, ct in enumerate(all_ct_list):
		this_ct_fn = all_ct_fn_list[ct_index] # get the overlap_enrichment file associated with this cell type
		this_ct_enrichment_df = get_one_enrichment_ct_model_df(this_ct_fn, num_state_ct_model) # get the enrichment data frame associated with this cell type
		max_enrich_value_all_ct_df[ct] = this_ct_enrichment_df['max_enrichment'] # store the values of maximum enrichment in this cell type
		max_enrich_25_state_all_ct_df[ct] = this_ct_enrichment_df['max_enrichment_context'] # store the names of the state that is most enriched with each of the full-stack state for this one cell type
	print ("Done getting data from all cell types " )
	get_rank_25_state_df(max_enrich_25_state_all_ct_df, max_enrich_value_all_ct_df, output_fn)
	print ("Done ranking enrichment data across cell types ")

get_ranked_max_enriched_25_state(args.input_folder, args.output_fn, args.num_state_ct_model)

