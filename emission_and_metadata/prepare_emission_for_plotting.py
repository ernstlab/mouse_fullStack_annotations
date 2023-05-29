import os
import pandas as pd 
import numpy as np 
emission_fn = './data/emissions_100.txt'
emission_df = pd.read_csv(emission_fn, header = 0, index_col = None, sep = '\t')
emission_df = emission_df.T # rows: experiments, columns: states
emission_df.reset_index(inplace =  True) # experiments are now a columns
emission_df.columns = emission_df.iloc[0] # first row is the states --> become headers
emission_df = emission_df.rename(columns={'State (Emission order)': 'experiments'})
emission_df = emission_df.drop(0, axis = 0) # after first row becomes headers, drop first row
# experiments are of form hindbrain_ATAC-seq_ENCSR662KNY
emission_df['assay'] = emission_df['experiments'].apply(lambda x: x.split('_')[-2].split('-')[0])
emission_df['biosample'] = emission_df['experiments'].apply(lambda x: x.split('_')[0])
emission_df['expID'] = emission_df['experiments'].apply(lambda x: x.split('_')[2])
assay_meta_fn = './data/assay_count.csv'
biosample_meta_fn = './data/biosample_term_name_counts.csv'
assay_df = pd.read_csv(assay_meta_fn, header = 0, index_col = None, sep = ',')[['mark', 'color', 'big_group']]
assay_df = assay_df.rename(columns = {'color': 'mark_color', 'big_group': 'assay_big_group'})
print(dict(assay_df.groupby(['mark', 'mark_color']).groups.keys()))
biosample_df = pd.read_csv(biosample_meta_fn, header = 0, index_col = None, sep = ',')[['biosample_name', 'group', 'big_group', 'color']]
biosample_df['group'] = biosample_df['group'].apply(lambda x: x.lower())
biosample_df = biosample_df.rename(columns = {'color': 'biosample_color', 'group': 'biosample_group', 'big_group': 'biosample_big_group', 'color': 'biosample_color'})
biosample_df['biosample_group'] = biosample_df['biosample_group'].replace('', 'other')
print(dict(biosample_df.groupby(['biosample_group', 'biosample_color']).groups.keys()))

meta_fn = './data/slims_metadata.txt'
meta_df = pd.read_csv(meta_fn, header = 0, index_col = None, sep = '\t') # 'experimentID', 'organ_slims', 'cell_slims', 'developmental_slims', 'system_slims', 'biosample_summary', 'simple_biosample_summary'
meta_df['organ_group'] = meta_df['organ_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('musculature_of_body', 'musculature').replace('', 'unknown')
meta_df['cell_group'] = meta_df['cell_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('', 'unknown')
meta_df['development_group'] = meta_df['developmental_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('', 'unknown') # The ectoderm gives rise to the skin and the nervous system. The mesoderm specifies the development of several cell types such as bone, muscle, and connective tissue. Cells in the endoderm layer become the linings of the digestive and respiratory system, and form organs such as the liver and pancreas
meta_df['system_group'] = meta_df['system_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('', 'unknown').apply(lambda x: x.split('_system')[0])

# PROCEESSING METADATA ABOUT THE ORGANS AND DEVELOPMENTAL STAGES
organ_meta_fn = './data/organ.txt'
organ_meta_df = pd.read_csv(organ_meta_fn, header = None, index_col = None, sep = '\t') # organ, color, order
ORGAN_GROUP_COLOR = dict(zip(organ_meta_df[0], organ_meta_df[1])) # keys: organ, values: color
print(ORGAN_GROUP_COLOR)
ORGAN_ORDER = dict(zip(organ_meta_df[0], organ_meta_df[2])) # keys: organ, values: organ numbered order
meta_df['organ_color'] = meta_df['organ_group'].apply(lambda x: ORGAN_GROUP_COLOR[x])
meta_df['organ_order'] = meta_df['organ_group'].apply(lambda x: ORGAN_ORDER[x])
meta_df.drop(labels = ['cell_slims', 'organ_slims', 'system_slims'], axis = 1, inplace = True)

print(emission_df.columns)
print(biosample_df.columns)
emission_df = emission_df.merge(assay_df, how = 'left', left_on = 'assay', right_on = 'mark')
emission_df = emission_df.merge(biosample_df, how = 'left', left_on = 'biosample', right_on = 'biosample_name')
print(emission_df[emission_df['biosample_group'].isnull()])
print
emission_df = emission_df.merge(meta_df, how = 'left', left_on = 'expID', right_on = 'experimentID')
emission_df['biosample_group'] = emission_df['biosample_group'].replace('', 'other')
save_fn = './data/emissions_100_for_pheatmap.txt'
emission_df.to_csv(save_fn, header = True, index = False, sep = '\t')
