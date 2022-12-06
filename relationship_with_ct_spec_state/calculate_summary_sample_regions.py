# this file will take output of code: count_sample_region_per_full_stack_state.py and produce some summary tables of the sample regions with the different seeds. 
# one summary table will have full_stack state as rows, each column will correspond to <cell_group>_<stateCT>, where stateCT is each of the 25 state. The table will be painted in excel and the color scale will correspond to the the color of the 25 states. It will be a beautiful excel. 
import os 
import sys
import pandas as pd 
import helper
import numpy as np 
import seaborn as sns
import glob
from scipy.stats import mannwhitneyu
import argparse
parser = argparse.ArgumentParser(description = 'Create an excel of the top ct-spec states enriched in each full-stack state')
parser.add_argument('--all_seed_folder', type=str,
	help = "where there are multiple subfolders, each containing enrichment data for different cell type specific model")
parser.add_argument('--output_folder', type=str,
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
helper.check_dir_exist(args.all_seed_folder)
helper.make_dir(args.output_folder)
helper.check_file_exist(args.full_state_annot_fn)
helper.check_file_exist(args.ct_state_annot_fn)
helper.check_file_exist(args.ct_group_fn)

NUM_FULL_STACK_state = 100
NUM_NON_CT_COLUMNS = 4 # the first NUM_NON_CT_COLUMNS in bed_map_df are chrom, start, end, full_stack

def get_ct_spec_state_annot(ct_state_annot_fn):
	df = pd.read_csv(ct_state_annot_fn, header = 0, index_col = None, sep = '\t') # state, mneumonics, meaning, color
	CT_STATE_MNEUNOMICS_DICT = pd.Series(df.mneumonics.values, index = df.state).to_dict() # keys: state(1-based), values: state mneumonics
	return CT_STATE_MNEUNOMICS_DICT, df

CT_STATE_MNEUNOMICS_DICT, CT_STATE_DF = get_ct_spec_state_annot(args.ct_state_annot_fn)

def get_full_stack_annot_df(full_state_annot_fn):
	fS_annot_df = pd.read_csv(full_state_annot_fn, header = 0, index_col = None, sep = '\t') # state, mneumonics, Long annotations, Short Annotations, group, color, state_order_by_group. state: 1--> 100, color: hex
	fS_annot_df = fS_annot_df[['state', 'mneumonics', 'color', 'state_order_by_group']]
	return fS_annot_df

def read_ct_annot_df(ct_group_fn):
	unique_cell_groups = ['forebrain', 'midbrain', 'hindbrain', 'neural-tube', 'facial-prominence', 'limb', 'intestine', 'stomach', 'liver', 'kidney', 'lung', 'heart']#CUSTOM: change this line of code for better customization
	unique_stages = ['e10.5', 'e11.5', 'e12.5', 'e13.5', 'e14.5', 'e15.5', 'e16.5', 'P0']#CUSTOM: change this line of code for better customization
	ct_annot_df = pd.read_csv(ct_group_fn, header = 0, index_col = None, sep = '\t')
	ct_annot_df = ct_annot_df.rename(columns = {'ct': 'cell_group'})
	cell_group_color_dict = dict(zip(ct_annot_df.cell_group, ct_annot_df.color))
	ct_annot_df = ct_annot_df[['tissue_stage', 'cell_group', 'stage']] #, 'COLOR', 'ANATOMY']] 
	ct_annot_df['cell_group'] = ct_annot_df.cell_group.astype('category')
	(ct_annot_df.cell_group).cat.set_categories(unique_cell_groups, inplace = True)
	ct_annot_df['stage'] = (ct_annot_df.stage).astype('category')
	(ct_annot_df.stage).cat.set_categories(unique_stages, inplace = True)
	ct_annot_df = ct_annot_df.sort_values(by = ['cell_group', 'stage']) # sort by gruups so that we can later arrange the data in specific order
	ct_annot_df['cell_group'] = ct_annot_df['cell_group'].astype(str)
	ct_annot_df['stage'] = ct_annot_df['stage'].astype(str)
	ct_annot_df = ct_annot_df.drop('stage', axis = 1)
	return ct_annot_df, cell_group_color_dict

def process_one_seedFile(bed_map_fn, num_state_ct_model): # we process the result of one file that contains the sampled reference 
	bed_map_df = pd.read_csv(bed_map_fn, header = 0, index_col = None, sep = '\t')
	# rearrange the columns such that cell types of the same group are together
	num_segment_per_state = bed_map_df.shape[0] / NUM_FULL_STACK_state
	# chrom, start, end, full_stack, E001, E002, etc. 
	input_ct_list = bed_map_df.columns[NUM_NON_CT_COLUMNS:] # the first 4 columns are chrom, start, end, full_stack, the following columns correspond to different cell types
	all_stateCT_list = list(map(lambda x: 'U' + str(x+1), range(num_state_ct_model))) # NOTE: this is a line of code that can be better automated such that the 'U' or 'E' is not hard codedE1 --> E25
	bed_map_df = bed_map_df.groupby('full_stack') 
	result_df = pd.DataFrame(columns = ['stateCT', 'ct', 'prop_in_ct', 'full_stack'])
	for full_stack_state, state_df in bed_map_df:# loop through each group of full-stack state
		state_df = state_df[input_ct_list] 
		count_df = state_df[input_ct_list].apply(pd.Series.value_counts, normalize = True) # only get columns corresponding to the cell types in ROADMAP, then count the number of occurrences of each of the 25 states in each ROADMPA cell type. Then, calculate the fraction of occurrence for each of the 25 state in each of the cell type (normalize = True). There may not be all 25 states, and there will be NAN values (later replace by 0). Rows: each of the 25 states. columns: ROADMAP cell types
		count_df = count_df.fillna(0) # Rows: each of the 25 states (may not be complete all 25 states). columns: ROADMAP cell types
		missing_state_list = np.setdiff1d(all_stateCT_list, count_df.index)
		count_df = count_df.append(pd.DataFrame(0, columns = count_df.columns, index = missing_state_list)) # add the missing states from the 25 states, and say that the proportion of these states in the sampled regions in each ct is 0
		count_df = count_df.reset_index().rename(columns = {'index' : 'stateCT'}) # one more column showing the state E1 --> E25
		count_df = count_df.melt(id_vars = ['stateCT']).rename(columns = {'variable' : 'ct', 'value' : 'prop_in_ct'}) # 3 columns: stateCT, ct, prop_in_ct showing the proportion that each stateCT is sampled in each ct from ROADMAP
		count_df['full_stack'] = full_stack_state # 4 columns: stateCT, ct, prop_in_ct , full_stack states
		result_df = result_df.append(count_df)
	return result_df

def calculate_avg_proportion_across_ct_per_group(oneSeed_prop_df, ct_annot_df):
	# ct_annot_df has 2 columns: tissue_stage and group
	# oneSeed_prop_df is a df with 4 columns stateCT, ct, prop_in_ct, full_stack
	oneSeed_prop_df = oneSeed_prop_df.merge(ct_annot_df, how = 'left', left_on = 'ct', right_on = 'tissue_stage') # 6 columns stateCT, ct, prop_in_ct, full_stack, tissue_stage, group
	oneSeed_prop_df = oneSeed_prop_df[oneSeed_prop_df['stateCT'] != '.'] # filter out lines with '.'
	piv_df = pd.pivot_table(oneSeed_prop_df, values = 'prop_in_ct', index = ['full_stack', 'stateCT'], columns = ['cell_group'], aggfunc = {'prop_in_ct' : np.mean}) # a pivot table index has two layers: full_stack and stateCT. Columns correspond to group
	piv_df = piv_df.unstack() # df has columns with two layers: group and stateCT, index: full_stack
	return piv_df

def color_full_stack_state(column_data, full_stack_color_dict):
	results = [''] * len(column_data)
	for row_index, value in enumerate(column_data):
		if row_index < NUM_FULL_STACK_state:
			try:	
				color = full_stack_color_dict[value]
			except:
				color = '#ffffff' # white
		else:
			color = '#ffffff' # white
		results[row_index] = 'background-color: %s' % color
	return results

def color_annot_cellgroup_stateCT(column_data, cell_group_color_dict, stateCT_color_dict): 
	results = [''] * len(column_data)
	for index, value in enumerate(column_data):
		if index == 0:
			try:
				color = cell_group_color_dict[value]
			except:
				color = '#ffffff' # white
		elif index == 1:
			try:
				color = stateCT_color_dict[value]
			except:
				color = '#ffffff' # white
		else:
			color = '#ffffff' # white
		results[index] = 'background-color: %s' % color
	return results

def color_avg_allSeed_group_df(avg_allSeed_group_df, CT_STATE_DF, num_state_ct_model):
	# test_df has exactly the same columns as avg_allSeed_group_df: mneumonics, <cell_group>_E<stateCT>, state (1 ,..., 100), color
	full_stack_color_dict = dict(zip(avg_allSeed_group_df.mneumonics, avg_allSeed_group_df.color))
	stateCT_color_dict = dict(zip(CT_STATE_DF.mneumonics, CT_STATE_DF.color))
	stateCT_color_dict['15_Ns'] = '#3498DB' 
	colored_df = avg_allSeed_group_df.style.apply(lambda x: color_full_stack_state(x, full_stack_color_dict), axis = 0) # color full_stack states column
	# now onto coloring individual states
	NUM_ROWS_TO_COLOR_stateCT_WISE = NUM_FULL_STACK_state + 1  # full stack states and min and max, stateCT-block wise
	NUM_ROWS_TO_COLOR_THEOR_WISE = NUM_FULL_STACK_state + 3 # full stack states and theoretical min (0) and max (1)
	for state_index in range(num_state_ct_model):
		state_str = 'U' + str(state_index + 1) # 0 --> E1, NOTE this can be replaced by something such that the U can be user-specified
		state_color = stateCT_color_dict[CT_STATE_DF.loc[state_index, 'mneumonics']]
		cm = sns.light_palette(state_color, as_cmap=True) # the color map for this stateCT
		colnames_to_colors = pd.Series(avg_allSeed_group_df.columns).str.endswith('_' + state_str) # a list of true or false. True if the column name ends with the state that we want to color
		colnames_to_colors = avg_allSeed_group_df.columns[colnames_to_colors] #now just a list of column names to pass into pd.IndexSlice
		for colname in colnames_to_colors:
			colored_df = colored_df.background_gradient(subset = pd.IndexSlice[:NUM_ROWS_TO_COLOR_stateCT_WISE, [colname]], cmap = cm)
	return colored_df

def helper_get_stateCT_from_colname(colname, CT_STATE_MNEUNOMICS_DICT):
	s_list = colname.split('_')
	if len(s_list) < 2:
		return ''
	else: 
		return CT_STATE_MNEUNOMICS_DICT[int(s_list[1][1:])]

def color_column_annotations(avg_df_columns, cell_group_color_dict, CT_STATE_DF): #avg_df_columns should be the columns of avg_allSeed_group_df
	stateCT_color_dict = dict(zip(CT_STATE_DF.mneumonics, CT_STATE_DF.color))
	stateCT_color_dict['15_Ns'] = '#3498DB'	 
	df = pd.DataFrame(columns = avg_df_columns)
	df.loc[0] = list(map(lambda x: x.split('_')[0], avg_df_columns))
	df.loc[1] = list(map(lambda x: helper_get_stateCT_from_colname(x, CT_STATE_MNEUNOMICS_DICT), avg_df_columns))
	colored_df = df.style.apply(lambda x: color_annot_cellgroup_stateCT(x, cell_group_color_dict, stateCT_color_dict), axis = 0)
	return colored_df

def calculate_stateCT_max_min_row(avg_allSeed_group_df, num_state_ct_model):
	"""
	For each stateCT, there are <#cell_group> columns, each column has 100 rows corrsponding to 100 full-stack state. In this function, we will report the maximum values of proportion in across a subset of avg_allSeed_group_df with nrow = # full_stack_state and ncols  = # cell_group, corresponding to each state. The result  is a list of length <#cell_group> * 25, where every <#cell_group> entry has the same value of maximum/minimum proportion
	"""
	results_max = []
	results_min = []
	for state_index in range(num_state_ct_model):
		state_str = 'U' + str(state_index + 1) # 0 --> U1 # NOTE this line of code can be made such that the U or the E can be change to user-specified
		colnames_to_this_state = pd.Series(avg_allSeed_group_df.columns).str.endswith('_' + state_str) # a list of true or false. True if the column name ends with the state that we want to color
		colnames_to_this_state = avg_allSeed_group_df.columns[colnames_to_this_state] #now just a list of column names to pass into pd.IndexSlice
		this_group_max = (avg_allSeed_group_df[colnames_to_this_state]).values.max(axis = None) # max value of the entire dataframe (this subset of the big dataframe)
		results_max += [this_group_max] * len(colnames_to_this_state)
		this_group_min = (avg_allSeed_group_df[colnames_to_this_state]).values.min(axis = None) # min value of the entire dataframe (this subset of the big dataframe)
		results_min += [this_group_min] * len(colnames_to_this_state)
	return results_max, results_min

def test_wilcoxon_group_spec(oneSeed_proportion_df_list, ct_annot_df, fS_annot_df, output_folder, num_state_ct_model):
	"""
	oneSeed_proportion_df_list from count_sample_region_each_full_stack_state: each item in the list is a df with 4 columns stateCT, ct, prop_in_ct, full_stack. Each df corresponds to the results from one random generator seed
	ct_annot_df has columns: tissue_stage, group
	Note: both full_stack and stateCT are now in the form E<state_index_1_based>
	"""
	allSeed_prop_df = pd.concat(oneSeed_proportion_df_list) # 1 big df with 4 columns stateCT, ct, prop_in_ct, full_stack
	# allSeed_prop_df = pd.read_csv('./trial.gz', header = 0, index_col = None, sep = '\t') # this is for debug reasons
	allSeed_prop_df = allSeed_prop_df.merge(ct_annot_df, how = 'left', left_on = 'ct', right_on = 'tissue_stage')
	allSeed_prop_df = allSeed_prop_df.drop(columns = ['ct', 'tissue_stage']) # drop the unnecessary columns
	uniq_group_list = ['forebrain', 'midbrain', 'hindbrain', 'neural-tube', 'facial-prominence', 'limb', 'intestine', 'stomach', 'liver', 'kidney', 'lung', 'heart']#CUSTOM: change this line of code for better customization
	output_columns = []
	for state in range(num_state_ct_model):
		output_columns += list(map(lambda x: x + "_E" + str(state+1), uniq_group_list)) 
	result_df = pd.DataFrame(columns = ['full_stack'] + output_columns) # we add full_stack but not include in output_columns because we want to use output_columns later when merge with fS_state_annot_df
	allSeed_prop_df = allSeed_prop_df.groupby('full_stack')
	for fS_state in range(NUM_FULL_STACK_state):
		state_str = 'E' + str(fS_state + 1)
		this_fS_df = allSeed_prop_df.get_group(state_str).groupby('stateCT') # this df contains all data associated with this full_stack state
		result_row = [fS_state+1] # full_stack_state, 1-based, int
		for stateCT in range(num_state_ct_model): # 25
			stateCT_str = 'U' + str(stateCT+1) # NOte this code can be rewritten to make the U more of an user choice
			this_s25_df = this_fS_df.get_group(stateCT_str) # stateCT: current stateCT, ct, prop_in_ct, full_stack: current fS state, group
			for cell_group in uniq_group_list:
				row_bool_in_cgroup = pd.Series(this_s25_df['cell_group']).str.match(cell_group)
				x = np.array(this_s25_df['prop_in_ct'][row_bool_in_cgroup]) # in this group
				y = np.array(this_s25_df['prop_in_ct'][~row_bool_in_cgroup]) # NOT in the group
				try:  
					t = mannwhitneyu(x, y , use_continuity = False, alternative = 'greater')
					result_row.append(t.pvalue)
				except: # if the test crash, that's because all values in x and y are identical --> p value should be 1
					result_row.append(1)
		result_df.loc[result_df.shape[0]] = result_row
	result_df = result_df.merge(fS_annot_df, how = 'left', left_on = 'full_stack', right_on = 'state')
	result_df = result_df.sort_values('state_order_by_group') 
	result_df.reset_index(drop = True, inplace = True) # reset index after reordering the rows 
	result_df = result_df[['mneumonics'] + output_columns + ['state', 'color']] # rearrange the columns to be similar to avg_allSeed_group_df
	save_fn = os.path.join(output_folder, 'pval_mannU_test_across_ct.csv.gz')
	result_df.to_csv(save_fn, header = True, index = False, sep = '\t', compression = 'gzip') # debug
	return result_df

def color_significant_pvalue(column_data, pval_thres):
	results = [''] * len(column_data)
	for row_index, value in enumerate(column_data):
		if isinstance(value, float) and value <= pval_thres:
			results[row_index] = 'background-color: #A7CCE5' # blue sky
	return results

def color_pvalue_df(test_df, pval_thres):
	print("pval_thres: " + str(pval_thres))
	full_stack_color_dict = dict(zip(test_df.mneumonics, test_df.color))
	colored_df = test_df.style.apply(lambda x: color_full_stack_state(x, full_stack_color_dict), axis = 0) # color full_stack states column
	colored_df = colored_df.apply(lambda x: color_significant_pvalue(x, pval_thres), axis = 0) # color each column
	return colored_df

def report_only_significant_cells(test_df, pval_stringent_thres = 1.0E-100):
	# test_df: mneumonics, <cell_group>_<stateCT:Esomething>, state, color
	# pval_stringent_thres is actually much samller than the bonferroni-corrected p-value threshold
	pval_df = test_df.iloc[:,1:-2] # skip the first and the last two columns beacuse those correspond to mneumonics, state and color
	rows_pick = pval_df.le(pval_stringent_thres).any(axis = 1) # T/F for each rows that has >=1 pvalues <= pval_stringent_thres
	columns_pick = pval_df.columns[pval_df.le(pval_stringent_thres).any(axis = 0)] # list of column names for columns that has >= 1 pvalue < pval_stringent_thres	
	uniq_group_list = np.unique(list(map(lambda x: x.split('_')[0], columns_pick)))
	rearranged_colnames = [] # rearrange such that columns of the same cell groups are put next to each other
	columns_pick = pd.Series(columns_pick) # transform to pandas Series to use some useful functions
	for group in uniq_group_list:
		colnames_this_group = columns_pick[columns_pick.str.startswith(group)]
		rearranged_colnames += list(colnames_this_group)
	significant_df = test_df.loc[rows_pick, ['mneumonics'] + list(rearranged_colnames) + ['state', 'color']] # select only rows and columns with the significant p values, add the columns corresponding to mneumonics, state, color
	return significant_df

def create_excel_from_precalculated_data(output_folder, all_seed_folder, full_state_annot_fn, ct_group_fn, num_state_ct_model):
	# This function was designed out of necessity for debuging purposes. We need this function after we have calculated the data for files pval_mannU_test_across_ct.csv.gz and avg_prop_stateCT_per_group.csv.gz through the function count_sample_region_each_full_stack_state, but the formatting kept being modified and we don't want to keep recalculating these files through the function count_sample_region_each_full_stack_state. So, to save time, I made this function. 
	ct_annot_df, cell_group_color_dict = read_ct_annot_df(ct_group_fn) # tissue_stage, group. We need the cell_group_color_dict for later use of coloring the excel file
	fS_annot_df = get_full_stack_annot_df(full_state_annot_fn) # 'state', 'mneumonics', 'color', 'state_order_by_group'
	test_df = pd.read_csv(os.path.join(output_folder, 'pval_mannU_test_across_ct.csv.gz'), header = 0, index_col = None, sep = '\t') 
	avg_allSeed_group_df = pd.read_csv(os.path.join(output_folder, "avg_prop_stateCT_per_group.csv.gz"), header = 0, index_col = None, sep = '\t') # this is for debug
	uniq_group_list = np.unique(list(map(lambda x: x.split('_')[0],avg_allSeed_group_df.columns[1:-2]))) # the first column is mneumonics, the last two columns are state and color 
	rearranged_colnames = []
	for stateCT in range(num_state_ct_model):
		rearranged_colnames += list(map(lambda x: x + '_U' + str(stateCT+1), uniq_group_list))
	# record the raw data, before we add some rows about the max and min values of each column
	num_columns_with_prop_values = len(rearranged_colnames) # number of columns that are associated with the numbers showing proportion of the sampled regions that are in the stateCT for each cell group 
	num_columns_with_empty_values = avg_allSeed_group_df.shape[1] - (1 + num_columns_with_prop_values) # number of columns that we will fill with empty at the end when we try to append some few rows before making the colored excel file. the 1 in the code corresponds to column 'full_stack'
	stateCT_wise_max, stateCT_wise_min = calculate_stateCT_max_min_row(avg_allSeed_group_df, num_state_ct_model)
	avg_allSeed_group_df.loc[avg_allSeed_group_df.shape[0]] = ['stateCT_wise_min'] + stateCT_wise_min + [''] * num_columns_with_empty_values # add one more row showing the min value in each block of stateCT in the table
	avg_allSeed_group_df.loc[avg_allSeed_group_df.shape[0]] = ['stateCT_wise_max'] + stateCT_wise_max + [''] * num_columns_with_empty_values # add one more row showing the max value in each block of stateCT in the table
	avg_allSeed_group_df.loc[avg_allSeed_group_df.shape[0]] = ['theoretical_min'] + [0] * num_columns_with_prop_values + [''] * num_columns_with_empty_values # add one more row showing the min value/ This is needed because we want to control the color of the cells in excel to be uniformly ranged (0,1) later
	avg_allSeed_group_df.loc[avg_allSeed_group_df.shape[0]] = ['theoretical_max'] + [1] * num_columns_with_prop_values + [''] * num_columns_with_empty_values # add one more row showing the max value/ This is needed because we want to control the color of the cells in excel to be uniformly ranged (0,1) later
	print ("Done getting all the necessary data for the excel")
	main_colored_df = color_avg_allSeed_group_df(avg_allSeed_group_df, CT_STATE_DF, num_state_ct_model)
	all_column_annot_colored_df = color_column_annotations(avg_allSeed_group_df.columns, cell_group_color_dict, CT_STATE_DF)
	ALPHA = 0.01
	bonreffroni_pval_thres = ALPHA / (num_state_ct_model * NUM_FULL_STACK_state * len(uniq_group_list)) # bonreffroni corrected p-value threshold
	pval_colored_df = color_pvalue_df(test_df, bonreffroni_pval_thres)
	pval_stringent_thres = 1.0E-100
	stringent_pval_df = report_only_significant_cells(test_df, pval_stringent_thres)
	significant_colored_df = color_pvalue_df(stringent_pval_df, pval_stringent_thres)
	significant_column_annot_df = color_column_annotations(stringent_pval_df.columns, cell_group_color_dict, CT_STATE_DF)
	output_fn = os.path.join(output_folder,  "avg_prop_stateCT_per_group.xlsx")
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	main_colored_df.to_excel(writer, sheet_name='avg_prop_stateCT_per_group')
	all_column_annot_colored_df.to_excel(writer, sheet_name = 'all_column_annot')
	pval_colored_df.to_excel(writer, sheet_name = 'pval_mwu_1s')
	significant_colored_df.to_excel(writer, sheet_name = 'pval_thres_1.0E-100')
	significant_column_annot_df.to_excel(writer, sheet_name = 'significant_column_annot')
	writer.save()
	print ("Done painting the excel")
	return

def count_sample_region_each_full_stack_state(all_seed_folder, output_folder, num_state_ct_model, full_state_annot_fn, ct_group_fn):
	ct_annot_df, cell_group_color_dict = read_ct_annot_df(ct_group_fn) # tissue_stage, group. We need the cell_group_color_dict for later use of coloring the excel file
	fS_annot_df = get_full_stack_annot_df(full_state_annot_fn) # 'state', 'mneumonics', 'color', 'state_order_by_group'
	bedmap_fn_list = glob.glob(all_seed_folder + 'seed_*/sample_segment_fullStack_ctState.bed.gz')
	oneSeed_proportion_df_list = list(map(lambda x: process_one_seedFile(x, num_state_ct_model), bedmap_fn_list)) # each item in the list is a df with 4 columns stateCT, ct, prop_in_ct, full_stack
	allSeed_prop_df = pd.concat(oneSeed_proportion_df_list) #for debug
	allSeed_prop_df.to_csv('trial.gz', header = True, index = False, sep = '\t', compression = 'gzip') # for debug
	test_df = test_wilcoxon_group_spec(oneSeed_proportion_df_list, ct_annot_df, fS_annot_df, output_folder, num_state_ct_model)
	oneSeed_proportion_df_list = list(map(lambda x: calculate_avg_proportion_across_ct_per_group(x, ct_annot_df), oneSeed_proportion_df_list)) # list of df with columns with two layers: group and stateCT, index: full_stack
	avg_allSeed_group_df = pd.concat(oneSeed_proportion_df_list).groupby(level = 0).mean()
	uniq_group_list =['forebrain', 'midbrain', 'hindbrain', 'neural-tube', 'facial-prominence', 'limb', 'intestine', 'stomach', 'liver', 'kidney', 'lung', 'heart']#CUSTOM: change this line of code for better customization
	print(uniq_group_list)
	avg_allSeed_group_df.columns = ['_'.join(col) for col in avg_allSeed_group_df.columns.values] # now columns are just one layer of the form 
	# after the following line of code, there is one more column to this dataframe
	avg_allSeed_group_df.reset_index(drop = False, inplace = True) # full_stack, Adipose_E1  Blood & T-cell_E1  Brain_E1 --> Sm. Muscle_E25  Thymus_E25  iPSC_E25
	avg_allSeed_group_df['full_stack'] = avg_allSeed_group_df['full_stack'].apply(lambda x: int(x[1:])) # E1 --> 1
	avg_allSeed_group_df = avg_allSeed_group_df.merge(fS_annot_df, how =  'left', left_on = 'full_stack', right_on = 'state').fillna('')
	avg_allSeed_group_df = avg_allSeed_group_df.sort_values('state_order_by_group')
	avg_allSeed_group_df.reset_index(drop = True, inplace = True) # after the sorting step, the indices are not in order so we just reset it.
	avg_allSeed_group_df = avg_allSeed_group_df.drop(columns = ['full_stack', 'state_order_by_group'])
	rearranged_colnames = []
	for stateCT in range(num_state_ct_model):
		rearranged_colnames += list(map(lambda x: x + '_U' + str(stateCT+1), uniq_group_list)) # Note this is something that we can change so that the U can be something else
	print(avg_allSeed_group_df.head())
	avg_allSeed_group_df = avg_allSeed_group_df[['mneumonics'] + rearranged_colnames + ['state', 'color']] # rearrange columns
	# record the raw data, before we add some rows about the max and min values of each column
	output_csv_fn = os.path.join(output_folder, "avg_prop_stateCT_per_group.csv.gz")
	avg_allSeed_group_df.to_csv(output_csv_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	print ("Done painting the excel")
	return

count_sample_region_each_full_stack_state(args.all_seed_folder, args.output_folder, args.num_state_ct_model, args.full_state_annot_fn, args.ct_group_fn)	
create_excel_from_precalculated_data(args.output_folder, args.all_seed_folder, args.full_state_annot_fn, args.ct_group_fn, args.num_state_ct_model)
print("Done!")

