import pandas as pd 
import numpy as np 
import glob
ct_meta_fn = '/Users/vuthaiha/Desktop/window_hoff/full_stacked_mouse/from_jason_061921/emissions/biosample_term_name_counts.csv'
ct_df = pd.read_csv(ct_meta_fn, header = 0, index_col = None, sep = ',')
lecif_folder = '/Users/vuthaiha/Desktop/window_hoff/data/ENCODE/mouse/chromHMM_15/from_lecif'
lecif_fn_list = glob.glob(lecif_folder + '/*_15_segments.bed.gz')
sample_list = list(map(lambda x: x.split('/')[-1].split('_15_segments.bed.gz')[0], lecif_fn_list))
result_df = pd.DataFrame({'tissue_stage' : sample_list})
result_df['ct'] = result_df['tissue_stage'].apply(lambda x: x.split('_')[1])
result_df['stage'] = result_df['tissue_stage'].apply(lambda x: x.split('_')[0])
result_df = result_df.merge(ct_df, how = 'left', left_on = 'ct', right_on = 'biosample_name')
result_df = result_df[['tissue_stage', 'ct', 'stage', 'color', 'group']]
facial_prominence_index = (result_df['ct'] == 'facial-prominence')
result_df.loc[facial_prominence_index, 'color'] = '#124c78' # dark blue
neural_tube_index = (result_df['ct'] == 'neural-tube')
result_df.loc[neural_tube_index, 'color'] = '#FFD924' # bright yellow
save_fn = '/Users/vuthaiha/Desktop/window_hoff/data/ENCODE/mouse/chromHMM_15/tissue_annotation.csv'
result_df.to_csv(save_fn, header = True, index = False, sep = '\t')
