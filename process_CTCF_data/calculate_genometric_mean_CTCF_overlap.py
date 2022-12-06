import pandas as pd 
import numpy as np 
import os
from scipy.stats import gmean 
from scipy.stats import gstd
import sys
import seaborn as sns
sys.path.append('/Users/vuthaiha/Desktop/window_hoff/source/mm10_annotations/overlap')
import create_excel_painted_overlap_enrichment_LECIF as lecif_code

fn = '/Users/vuthaiha/Desktop/window_hoff/full_stacked_mouse/from_jason_061921/overlap/ctcf/overlap_ctcf_mm10.txt'
save_excel_fn = '/Users/vuthaiha/Desktop/window_hoff/full_stacked_mouse/from_jason_061921/overlap/ctcf/overlap_ctcf_mm10.xlsx'
save_csv_fn = '/Users/vuthaiha/Desktop/window_hoff/full_stacked_mouse/from_jason_061921/overlap/ctcf/overlap_ctcf_mm10_with_calculated_mean.txt'
state_annot_fn = '/Users/vuthaiha/Desktop/window_hoff/full_stacked_mouse/from_jason_061921/state_annotation_processed.csv'
df = pd.read_csv(fn, header = 0, index_col = None, sep = '\t')
percent_genome_of_cont = df.loc[df.shape[0]-1,:]
df.drop(df.tail(1).index,inplace=True)
enr_colname = df.columns[df.columns.str.endswith('.bed.gz')]
print((df[enr_colname] == 0).sum())
def calculate_gmean(x):
	x = x.astype(float)
	non_zero_min = np.min(x[x!=0])
	x[x==0] = non_zero_min
	return gmean(np.array(x))

def calculate_geostd(x):
	x = x.astype(float)
	non_zero_min = np.min(x[x!=0])
	x[x==0] = non_zero_min
	return gstd(np.array(x))

df['geomean'] = df.apply(lambda x: calculate_gmean(x[enr_colname]), axis = 1)
df['geostd'] = df.apply(lambda x: calculate_geostd(x[enr_colname]), axis = 1)
df['mean'] = df.apply(lambda x: np.mean(x[enr_colname].astype(float)), axis = 1)
df['std'] = df.apply(lambda x: np.std(x[enr_colname].astype(float)), axis = 1)

cm = sns.light_palette("red", as_cmap=True)
try:
    df = df.rename(columns = {u'state (Emission order)': 'state', u'Genome %' : 'percent_in_genome', u'State (Emission order)' : 'state'})
except:
    pass
df = df.fillna(0) # if there are nan enrichment values, due to the state not being present (such as when we create files with foreground and background), we can fill it by 0 so that the code to make colorful excel would not crash.
(num_state, num_enr_cont) = (df.shape[0], df.shape[1] - 1)
enr_colname = list(map(lecif_code.get_rid_of_stupid_file_tail, enr_colname)) # fix the name of the enrichment context. If it contains the tail '.bed.gz' then we get rid of it
df.columns = list(map(lecif_code.get_rid_of_stupid_file_tail, df.columns))
# now change the column names of the enrichment contexts. If the contexts' names contain '.bed.gz' then get rid of it.
state_annot_df = lecif_code.get_state_annot_df(state_annot_fn)
df['state'] = (df['state']).astype(str).astype(int)
df = df.merge(state_annot_df, how = 'left', left_on = 'state', right_on = 'state', suffixes = ('_x', '_y'))
df = df.sort_values(by = 'state_order_by_group')
df = df.reset_index(drop = True)
for enr_cont in enr_colname:
    df[enr_cont] = df[enr_cont].astype(str).astype(float)  
num_remaining_columns = len(df.columns) - 2 - len(enr_colname)
df.loc[df.shape[0]] = list(percent_genome_of_cont) + [None] * num_remaining_columns  
mneumonics_index_in_row = df.columns.get_loc('mneumonics') # column index, zero-based of the menumonics entries, which we will use to paint the right column later
colored_df = df.style.background_gradient(subset = pd.IndexSlice[:(num_state-1), [enr_colname[0]]], cmap = cm)
for enr_cont in enr_colname[1:] + ['geomean', 'mean']:
    colored_df = colored_df.background_gradient(subset = pd.IndexSlice[:(num_state-1), [enr_cont]], cmap = cm)
colored_df = colored_df.apply(lambda x: lecif_code.color_state_annotation(x, [mneumonics_index_in_row]), axis = 1)
colored_df.to_excel(save_excel_fn, engine = 'openpyxl')
df.to_csv(save_csv_fn, header = True, index = False, sep = '\t')
