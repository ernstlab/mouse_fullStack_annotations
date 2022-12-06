import pandas as pd 
import numpy as np 
import seaborn as sns
import sys
sys.path.append('/u/home/h/havu73/project-ernst/source/mm10_annotations/overlap')
import create_excel_painted_overlap_enrichment_LECIF as lecif_code
from scipy import stats

def transform_state_to_number(state):
	'''
	E1 --> 1
	'''
	if state.startswith('E'):
		return int(state[1:])
	return int(state)
# 
def get_enrichment_df(enrichment_fn):
	enrichment_df = pd.read_csv(enrichment_fn, sep = "\t", header = 0)
	enrichment_df = enrichment_df.rename(columns = {u'state (Emission order)': 'state', u'Genome %' : 'percent_in_genome', u'State (Emission order)' : 'state'})
	enrichment_df = enrichment_df.fillna(0) # if there are nan enrichment values, due to the state not being present (such as when we create files with foreground and background), we can fill it by 0 so that the code to make colorful excel would not crash.
	(num_state, num_enr_cont) = (enrichment_df.shape[0] - 1, enrichment_df.shape[1] - 1)
	enr_cont_list = enrichment_df.columns[2:]
	enr_cont_list = list(map(lecif_code.get_rid_of_stupid_file_tail, enr_cont_list)) # fix the name of the enrichment context. If it contains the tail '.bed.gz' then we get rid of it
	enrichment_df.columns = list(enrichment_df.columns[:2]) + list(enr_cont_list)
	percent_genome_of_cont = enrichment_df.iloc[num_state, 2:]
	enrichment_df = enrichment_df.loc[:(num_state-1)] # rid of the last row because that is the row about percentage each genomic context occupies    
	# now change the column names of the enrichment contexts. If the contexts' names contain '.bed.gz' then get rid of it.
	no_state_df = enrichment_df[enrichment_df.columns[2:]] # only get the data of enrichment contexts for now, don't consider state and percent_in_genome
	enrichment_df['max_fold_context'] = no_state_df.apply(lambda x: x.idxmax(), axis = 1) # name of the context that are most enriched in this state (only applied for lecif score context here)
	return enrichment_df

def get_excel_highlight_lecif_phastCons_states(lecif_df, save_fn):
	'''
	lecif_df has columns state, lecif, phastConsEle60w, and columns related to state annotations
	'''
	lecif_top_states =  list((lecif_df[lecif_df['lecif'] == 1.0])['state'])
	print(len(lecif_top_states))
	phastCons_top_states = list(lecif_df.nlargest(16, columns = 'phastConsEle60w')['state'])
	print(len(phastCons_top_states))
	highlight_states = np.union1d(lecif_top_states, phastCons_top_states)
	lecif_df = lecif_df[lecif_df['state'].isin(highlight_states)]
	lecif_df = lecif_df.sort_values(['lecif', 'phastConsEle60w', 'state_order_by_group'], ascending = [True, False, True])
	cm = sns.light_palette("red", as_cmap=True)
	mneumonics_index_in_row = lecif_df.columns.get_loc('mneumonics') # column index, zero-based of the menumonics entries, which we will use to paint the right column later
	comment_index_in_row = lecif_df.columns.get_loc('comments')
	colored_df = lecif_df.style.apply(lambda x: lecif_code.color_state_annotation(x, [mneumonics_index_in_row, comment_index_in_row]), axis = 1)
	colored_df = colored_df.applymap(lecif_code.color_lecif_score, subset = pd.IndexSlice[:, ['lecif']])
	colored_df = colored_df.background_gradient(subset = pd.IndexSlice[:, ['phastConsEle60w']], cmap = cm)
	colored_df.to_excel(save_fn, engine = 'openpyxl')
	return 

lecif_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/overlap/lecif/avg_lecif_scores.txt'
gc_fn = '~/project-ernst/full_stacked_mouse/from_jason_061921/overlap/overlap_genome_context.txt'
state_annot_fn = '~/project-ernst/full_stacked_mouse/from_jason_061921/state_annotation_processed.csv'
save_fn = '~/project-ernst/full_stacked_mouse/from_jason_061921/overlap/highlight_states_lecif_phastcons.xlsx'

lecif_df = pd.read_csv(lecif_fn, header = 0, index_col = None, sep = '\t')
lecif_df['state'] = lecif_df['state'].apply(transform_state_to_number)
lecif_df = (lecif_df[['state', 'mean']]).rename(columns = {'mean': 'lecif'})

gc_df = get_enrichment_df(gc_fn)
gc_df['state'] = gc_df['state'].astype(int)
gc_df = gc_df[['state', 'phastConsEle60w']]

state_annot_df = lecif_code.get_state_annot_df(state_annot_fn)

lecif_df = lecif_df.merge(gc_df, left_on = 'state', right_on = 'state')
lecif_df.loc[:, 'state'] = lecif_df['state'].astype(int)
lecif_df = lecif_df.merge(state_annot_df, left_on = 'state', right_on = 'state')

print("Pearson correlation between lecif score ranges and phastConsEle60w fold enrichment: ")
print(np.corrcoef(lecif_df.lecif, lecif_df.phastConsEle60w))
print(stats.pearsonr(lecif_df.lecif, lecif_df.phastConsEle60w))
print("Spearman correlation between lecif score ranges and phastConsEle60w fold enrichment: ")
print(stats.spearmanr(lecif_df.lecif, lecif_df.phastConsEle60w))
print('10 states with highest enrichments with phastConsEle60w: ')
print('lecif_1.0')
# lecif1_df = lecif_df[lecif_df['lecif'] == 1.0].sort_values('state_order_by_group').copy()
# lecif1_df.reset_index(drop = True, inplace = True)
# get_excel_highlight_lecif_phastCons_states(lecif_df, save_fn)