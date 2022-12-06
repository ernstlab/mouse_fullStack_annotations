import pandas as pd 
import numpy as np 
def get_encode_metadata():
	# this file will take ing metadata_fn from ENCODE and emission_fn from Jason to produce a metadata file with relevant information for the emission matrix
	metadata_fn = '/u/home/h/havu73/project-ernst/data/ENCODE/Chip_peaks/mm10/metadata.csv'
	meta_df = pd.read_csv(metadata_fn, header = 0, index_col = None, sep = '\t')
	return meta_df

def get_metadata_from_jason():
	jason_meta_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/emissions/metadata_from_jason.txt'
	jasonM_df = pd.read_csv(jason_meta_fn, header = 0, index_col = None, sep = '\t')
	jasonM_df['organ_group'] = jasonM_df['organ_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('musculature_of_body', 'musculature').replace('', 'unknown')
	jasonM_df['cell_group'] = jasonM_df['cell_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('', 'unknown')
	jasonM_df['development_group'] = jasonM_df['developmental_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('', 'unknown') # The ectoderm gives rise to the skin and the nervous system. The mesoderm specifies the development of several cell types such as bone, muscle, and connective tissue. Cells in the endoderm layer become the linings of the digestive and respiratory system, and form organs such as the liver and pancreas
	jasonM_df['system_group'] = jasonM_df['system_slims'].apply(lambda x: '_'.join(x[1:-1].split(',')[0][1:-1].split())).replace('', 'unknown').apply(lambda x: x.split('_system')[0])
	jasonM_df = jasonM_df[['experimentID', 'organ_group', 'cell_group', 'development_group', 'system_group']]
	return jasonM_df

def fix_experiment_target(row):
	# if the experiment target is empty it's usually DNase or ATAC, which is annoted in the Assay column
	if pd.isnull(row['Experiment target']):
		return row['Assay']
	else:
		return row['Experiment target']

#############################################################################################
if __name__ == '__main__':
	emission_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/emissions/emissions_100.txt'
	emission_df = pd.read_csv(emission_fn, header = 0, index_col = None, sep = '\t')
	exp_code_list = pd.Series(emission_df.columns[1:]).apply(lambda x: x.split('_')[-1])
	exp_df = pd.DataFrame({'exp_id': exp_code_list})
	meta_df = get_encode_metadata()
	print('meta_df shape: {}'.format(meta_df.shape))
	jasonM_df = get_metadata_from_jason()
	print('jasonM_df shape: {}'.format(jasonM_df.shape))
	exp_df = exp_df.merge(jasonM_df, how = 'left', left_on = 'exp_id', right_on = 'experimentID')
	exp_df = exp_df.merge(meta_df, how = 'left', left_on = 'exp_id', right_on = 'Experiment accession')
	print('exp_df shape: {}'.format(exp_df.shape))
	exp_df = exp_df[['exp_id', 'organ_group', 'cell_group', 'development_group', 'system_group', 'Biosample term id','Biosample term name', 'Biosample type', 'Biosample organism', 'Biosample treatments', 'Biosample treatments amount', 'Biosample treatments duration', 'Biosample genetic modifications methods', 'Biosample genetic modifications categories', 'Biosample genetic modifications targets', 'Biosample genetic modifications gene targets', 'Biosample genetic modifications site coordinates', 'Biosample genetic modifications zygosity', 'Assay', 'Experiment target', 'File assembly', 'Donor(s)']]
	exp_df['Experiment target'] = exp_df.apply(fix_experiment_target, axis = 1) # beacuse
	grouped_df = exp_df.groupby('exp_id')

	result_df = pd.DataFrame(columns = exp_df.columns)
	for gName, gdf in grouped_df:
		if (gdf.nunique(axis = 0, dropna = True) > 1).any():
			# the count of unique values in each column in each gdf should be at most 1
			print('exp_df {} has non-unique values of other identifiers'.format(gName))
			continue
		elif(gdf['Experiment target'].isnull()).all():
			# print('exp_df {} has nan experiment target'.format(gName)) # It is because this expeirment is ATAC-seq or DNase-seq
			pass
		else:
			df = gdf.reset_index(drop = True)
			result_df.loc[result_df.shape[0]] = df.loc[0]


	result_df.dropna(axis = 1, how = 'all', inplace = True)
	output_fn = '/u/home/h/havu73/project-ernst/full_stacked_mouse/from_jason_061921/emissions/metadata.csv'
	result_df.to_csv(output_fn, header = True, index = False, sep = ',')