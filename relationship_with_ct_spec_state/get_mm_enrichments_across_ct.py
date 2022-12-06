import pandas as pd
import seaborn as sns
import numpy as np
import os
import helper
import glob
import argparse
parser = argparse.ArgumentParser(description = 'Create an excel of the top ct-spec states enriched in each full-stack state')
parser.add_argument('--input_folder', type=str,
	help = "where there are multiple subfolders, each containing enrichment data for different cell type specific model")
parser.add_argument('--output_fn', type=str,
	help = 'Where the files showing max, mean, min fold enrichment across all ct-state states are reported for all cell types for each state is reported')
parser.add_argument('--num_state_ct_model', type=int,
	help = 'number of states in the ct-spec models')
parser.add_argument('--full_state_annot_fn', type=str,
	help = 'file of the full-stack state annotations')
parser.add_argument('--ct_state_annot_fn', type=str,
	help = 'file of the ct-spec state annotations')
args = parser.parse_args()
print(args)
helper.check_dir_exist(args.input_folder)
helper.create_folder_for_file(args.output_fn)
helper.check_file_exist(args.full_state_annot_fn)
helper.check_file_exist(args.ct_state_annot_fn)

def get_ct_spec_state_annot(ct_state_annot_fn):
	df = pd.read_csv(ct_state_annot_fn, header = 0, index_col = None, sep = '\t') # state, mneumonics, meaning, color
	CT_STATE_MNEUNOMICS_LIST = list(df.mneumonics)
	CT_STATE_COLOR_DICT = pd.Series(df.color.values, index = df.mneumonics).to_dict() # keys: state mneumonics, values: state color
	return CT_STATE_MNEUNOMICS_LIST, CT_STATE_COLOR_DICT

CT_STATE_MNEUNOMICS_LIST, CT_STATE_COLOR_DICT = get_ct_spec_state_annot(args.ct_state_annot_fn)

def get_full_state_annot_df(full_state_annot_fn):
	df = pd.read_csv(full_state_annot_fn, header = 0, index_col = None, sep = '\t')
	FULL_STATE_COLOR_DICT = pd.Series(df.color.values, index = df.mneumonics)
	df = df[['state', 'mneumonics', 'state_order_by_group']]
	return FULL_STATE_COLOR_DICT, df

FULL_STATE_COLOR_DICT, FULL_STATE_ANNOT_DF = get_full_state_annot_df(args.full_state_annot_fn) #FULL_STATE_ANNOT_DF: state, mneumonics, state_order_by_group, FULL_STATE_COLOR_DICT: keys mneumonics, values color

def color_full_state_names(val):
	if val == "":
		color = '#FFFFFF'
	else:
		color = FULL_STATE_COLOR_DICT[val]
	return 'background-color: %s' % color

def get_rid_of_stupid_file_tail(context_name):
	if context_name.endswith('.bed.gz'):
		return(context_name[:-7])
	else:
		return(context_name)


def get_one_enrichment_ct_model_df(fn):
	df = pd.read_csv(fn, sep = '\t', header = 0)
	df = df.rename(columns = {'state (Emission order)' : 'state', 'Genome %': 'percent_in_genome'}) # rename some columns so that it is easier to write column names later
	context_name_list = list(map(get_rid_of_stupid_file_tail, df.columns[2:])) # get the names of enrichment contexts. If the enrichment contexts is <context>.bed.gz --> <context>
	df.columns = ['state', 'percent_in_genome'] + context_name_list
	df = df.drop([df.shape[0] - 1], axis = 0) # drop that last row, which is about the percentage of genome that each enrichment context occupies
	return df

def get_mmm_enrichment_one_genome_context(all_ct_df_list, all_ct_list, state_context_colName): 
	# state_context_colName: name of the column that we are looking at so that we know what column to look at for each data frame
	this_context_df = pd.DataFrame({'state' : (all_ct_df_list[0])['state']})
	# collect data from each of the cell type specific model --> columsn: cell type specific model, rows: 100 full stack states
	for df_index, df in enumerate(all_ct_df_list):
		ct_name = all_ct_list[df_index] # all_ct_fn_list and all_ct_df_list are responsive to each other, i.e. each element in each list corresponds to the same cell type. We tried zipping teh two list, but that gave a bug message. so this code here is not the most graceful.
		this_context_df[ct_name] = df[state_context_colName]
	this_context_df['max_enrichment'] = (this_context_df.drop('state', axis = 1)).max(axis = 1)
	this_context_df['min_enrichment'] = (this_context_df.drop('state', axis = 1)).min(axis = 1)
	this_context_df['median_enrichment'] = (this_context_df.drop('state', axis = 1)).median(axis = 1)
	return this_context_df

def get_25_state_annot(state):
	# '1' --> ''1_TssA''
	state_index = int(state) - 1 # zero-based
	return CT_STATE_MNEUNOMICS_LIST[state_index]

def get_CT_STATE_COLOR_DICT(state):
	# '1_TssA' --> 'red'
	state_index = int(state.split('_')[0]) # one-based
	color = CT_STATE_COLOR_DICT[state_index]
	return 'background-color: %s' % color

def paint_excel_mmm_enrichment(enrichment_df, num_state_ct_model):
	# here enrichment_df is actually max_enrichment_df from get_all_ct_model_mmm_enrichment
	cm = sns.light_palette("red", as_cmap=True)
	(num_state, num_enr_cont) = (enrichment_df.shape[0], enrichment_df.shape[1] - 1)
	state_colName_list = list(map(lambda x: 'U' + str(x + 1), range(num_state_ct_model))) #CUSTOM: this is a line of code that could be improved for better customization
	enrichment_df['state'] = enrichment_df['state'].astype(int)
	enrichment_df = enrichment_df.merge(FULL_STATE_ANNOT_DF, how = 'left', left_on = 'state', right_on = 'state')
	enrichment_df = enrichment_df[['state', 'percent_in_genome'] + state_colName_list + ['mneumonics', 'state_order_by_group']] # get the data frame to display columns in the expected order
	enrichment_df.columns = ['state', 'percent_in_genome'] + CT_STATE_MNEUNOMICS_LIST + ['mneumonics', 'state_order_by_group'] # rename the columns so that instead of 'state1.bed.gz' we have '1_TssA'
	enrichment_df = enrichment_df.sort_values(by = 'state_order_by_group')
	colored_df = enrichment_df.style.background_gradient(subset = pd.IndexSlice[:, CT_STATE_MNEUNOMICS_LIST], cmap = cm) # color the enrichment values into a red-white gradient in the first enrichment contnext
	colored_df = colored_df.applymap(color_full_state_names, subset = pd.IndexSlice[:, ['mneumonics']])
	return colored_df


def get_all_ct_model_mmm_enrichment(input_folder, output_fn, num_state_ct_model):
	all_ct_fn_list = glob.glob(input_folder + '/*/overlap_enrichment.txt')
	all_ct_list = list(map(lambda x: x.split('/')[-2], all_ct_fn_list)) # from '/path/to/E129/overlap_enrichment.txt' to 'E129' 
	all_ct_df_list = list(map(lambda x: get_one_enrichment_ct_model_df(x), all_ct_fn_list))
	all_ct_df_list = list(all_ct_df_list)
	# up until here, we have finished getting all the data from all the overlap enrichment files of all cell types into data frame --> all_ct_df_list [index] --> data frame of the cell type
	print ("Done getting all enrichment data from all cell types")
	# create a data frame that will report the median enrichment over all the cell types
	max_enrichment_df = pd.DataFrame({'state': (all_ct_df_list[0])['state']}) # this will report for each state in the 25-state model, and for each enrichment context, what are the highest enrichment in each state-context combination. 
	# now loop through each of the context that we care about
	median_enrichment_df = pd.DataFrame({'state': (all_ct_df_list[0])['state']}) # this will report for each state in the 25-state model, and for each enrichment context, what are the highest enrichment in each state-context combination. 
	# now loop through each of the context that we care about
	min_enrichment_df = pd.DataFrame({'state': (all_ct_df_list[0])['state']}) # this will report for each state in the 25-state model, and for each enrichment context, what are the highest enrichment in each state-context combination. 
	# now loop through each of the context that we care about
	enrichment_context_list = ['percent_in_genome'] + list((all_ct_df_list[0]).columns[2:]) # this is because we also want to report the max, median, min of percent_in_genome_of data from all ct-spec enrichments
	print (enrichment_context_list)
	for enr_context in enrichment_context_list:
		this_context_mmm_enrichment_df = get_mmm_enrichment_one_genome_context(all_ct_df_list, all_ct_list, enr_context)
		max_enrichment_df[enr_context] = this_context_mmm_enrichment_df['max_enrichment']
		min_enrichment_df[enr_context] = this_context_mmm_enrichment_df['min_enrichment']
		median_enrichment_df[enr_context] = this_context_mmm_enrichment_df['median_enrichment']		
	print ("Done getting all the necessary data!")
	max_colored_df = paint_excel_mmm_enrichment(max_enrichment_df, num_state_ct_model)
	min_colored_df = paint_excel_mmm_enrichment(min_enrichment_df, num_state_ct_model)
	median_colored_df = paint_excel_mmm_enrichment(median_enrichment_df, num_state_ct_model)
	# now save into 3 sheets
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	max_colored_df.to_excel(writer, sheet_name='max')
	min_colored_df.to_excel(writer, sheet_name='min')
	median_colored_df.to_excel(writer, sheet_name='median')
	writer.save()
	print ("Done getting the mmm_enrichment_df!")

get_all_ct_model_mmm_enrichment(args.input_folder, args.output_fn, args.num_state_ct_model)
