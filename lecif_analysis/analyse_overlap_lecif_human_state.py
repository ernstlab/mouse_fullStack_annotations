import pandas as pd 
import sys
sys.path.append('/u/home/h/havu73/project-ernst/source/mm10_annotations/overlap')
import create_excel_painted_overlap_enrichment_LECIF as lecif_code
import seaborn as sns
from matplotlib import pyplot as plt
import upsetplot # https://github.com/jnothman/UpSetPlot

def transform_state_to_number(state):
	'''
	E1 --> 1
	'''
	if state.startswith('E'):
		return int(state[1:])
	return int(state)

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
	return enrichment_df


def get_lecif_df():
	lecif_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/overlap/lecif/avg_lecif_scores.txt'
	lecif_df = pd.read_csv(lecif_fn, header = 0, index_col = None, sep = '\t')
	lecif_df['state'] = lecif_df['state'].apply(transform_state_to_number)
	lecif_df = (lecif_df[['state', 'mean']]).rename(columns = {'mean': 'lecif'})
	return lecif_df

def get_phastCons_df():
	gc_fn = '~/project-ernst/full_stacked_mouse/from_jason_061921/overlap/overlap_genome_context.txt'
	gc_df = get_enrichment_df(gc_fn)
	gc_df['state']= gc_df['state'].astype(int)
	gc_df = gc_df[['state', 'phastConsEle60w']]
	return gc_df


def color_excel(df, save_fn):
	print(df.columns)
	red_cm = sns.light_palette("red", as_cmap=True)
	blue_cm = sns.light_palette('blue', as_cmap = True)
	colored_df = df.style.background_gradient(subset = pd.IndexSlice[:, 'max_human_enr'], cmap = red_cm)
	colored_df = colored_df.background_gradient(subset = pd.IndexSlice[:, 'phastConsEle60w'], cmap = red_cm)
	colored_df = colored_df.background_gradient(subset = pd.IndexSlice[:, 'lecif'], cmap = blue_cm)
	indices_to_color_state = [df.columns.get_loc(x) for x in ['mneumonics', 'all',  'lecif_human', 'lecif_phastC', 'human_phastC', 'lecif_only', 'human_only', 'phastC_only']] # column index, zero-based of the mneumonics entries, which we will use to paint the right column later
	colored_df  = colored_df.apply(lambda x: lecif_code.color_state_annotation(x, indices_to_color_state), axis = 1)
	colored_df.to_excel(save_fn, engine = 'openpyxl')
	return

def report_state_satisfy_rank(row, include_rank_list, exclude_rank_list, n):
	'''
	return the state if the ranks in include_rank_list is <= n AND the ranks in exclude_rank_list is > n
	'''
	for rank_colname in exclude_rank_list:
		if row[rank_colname] <= n:
			return ''
	for rank_colname in include_rank_list:
		if row[rank_colname] > n:
			return ''
	return row['mneumonics']

def get_state_lists_based_on_rank(df, rank_colname, n):
	filter_df = df[df[rank_colname]  <= n]
	return list(filter_df['mneumonics'])

def get_upset_plot(df, save_fn, n):
	human_state = get_state_lists_based_on_rank(df, 'rank_human', n)
	lecif_state = get_state_lists_based_on_rank(df, 'rank_lecif', n)
	phastCons_state  = get_state_lists_based_on_rank(df, 'rank_phastCons', n)
	upset_content = upsetplot.from_contents({'human_state': human_state, 'LECIF': lecif_state, 'phastConsEle60w': phastCons_state})
	upsetplot.UpSet(upset_content, subset_size='count', show_counts='%d').plot()
	# print(plt)
	plt.savefig(save_fn, format='pdf', dpi=None, bbox_inches='tight')
	return

if __name__ == '__main__':
	mouse_state_annot_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/state_annotation_processed.csv'
	mouse_sAnnot_df = lecif_code.get_state_annot_df(mouse_state_annot_fn)
	human_fullStack_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/overlap/universalmouse_humanmapping.csv'
	df = pd.read_csv(human_fullStack_fn, header = 0, index_col = None, sep = ',')
	df.dropna(axis = 1, how = 'all', inplace = True)
	df = df.rename(columns= {'Mouse State # (Emission ordeing)': 'state', 'Most Enriched Human State' :'human_state', 'Enrichment': 'max_human_enr'})
	lecif_df = get_lecif_df()
	df = df.merge(lecif_df, left_on = 'state', right_on = 'state')
	df = df[['state', 'human_state', 'max_human_enr', 'lecif']]
	phastCons_df = get_phastCons_df()
	df = df.merge(phastCons_df, left_on = 'state', right_on = 'state')
	df = df.merge(mouse_sAnnot_df, left_on = 'state', right_on = 'state')
	df['rank_lecif'] = df['lecif'].rank(ascending = False).astype(int)
	df['rank_human'] = df['max_human_enr'].rank(ascending = False).astype(int)
	df['rank_phastCons'] = df['phastConsEle60w'].rank(ascending = False).astype(int)
	n = 20
	df['all'] = df.apply(lambda x: report_state_satisfy_rank(x, ['rank_lecif', 'rank_human', 'rank_phastCons'], [], n), axis = 1)
	df['lecif_human'] = df.apply(lambda x: report_state_satisfy_rank(x, ['rank_lecif', 'rank_human'], ['rank_phastCons'], n), axis = 1)
	df['lecif_phastC'] = df.apply(lambda x: report_state_satisfy_rank(x, ['rank_lecif', 'rank_phastCons'], ['rank_human'], n), axis = 1)
	df['human_phastC'] = df.apply(lambda x: report_state_satisfy_rank(x, ['rank_human', 'rank_phastCons'], ['rank_lecif'], n), axis = 1)
	df['lecif_only'] = df.apply(lambda x: report_state_satisfy_rank(x, ['rank_lecif'], ['rank_human', 'rank_phastCons'], n), axis = 1)
	df['human_only'] = df.apply(lambda x: report_state_satisfy_rank(x, ['rank_human'], ['rank_lecif', 'rank_phastCons'], n), axis = 1)
	df['phastC_only'] = df.apply(lambda x: report_state_satisfy_rank(x, ['rank_phastCons'],  ['rank_lecif', 'rank_human'], n), axis = 1)
	df = df.sort_values('state_order_by_group')
	excel_save_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/overlap/lecif/avg_lecif_phastCons_human.xlsx'
	color_excel(df, excel_save_fn)
	upset_save_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/overlap/lecif/upset_lecif_human_phastCons.pdf'
	get_upset_plot(df, upset_save_fn, n)
