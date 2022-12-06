import pandas as pd 
import numpy as np 
import helper 
import os 
import argparse
import sys
import requests
parser = argparse.ArgumentParser(description = 'This file will filter out the data needed to get the CTCF peaks in mouse')
parser.add_argument('--metadata_fn', type = str, required = True, help = 'metadata file that I got from ENCODE: mouse, CTCF, TF Chip-seq')
parser.add_argument('--download', action='store_true', default = False, help = 'If this flag is present, we will download the data of CTCF peaks')
parser.add_argument('--download_folder', type = str, required = True, help = 'Where data should be downloaded to')
parser.add_argument('--output_fn', required = True, help = 'The file of metadata that we will show as additional data for the paper')
args = parser.parse_args()
helper.check_file_exist(args.metadata_fn)
helper.create_folder_for_file(args.output_fn)
helper.make_dir(args.download_folder)
###################################################################################
def cleanup_metadata(metadata_fn, output_fn):
	meta_df = pd.read_csv(metadata_fn, header = 0, index_col = None, sep = '\t')
	meta_df.dropna(axis = 1, how = 'all',inplace = True) # drop all empty columns, only 34 of them left
	meta_df = meta_df[meta_df['File assembly'] == 'mm10']
	meta_df = meta_df[meta_df['Assay'] == 'TF ChIP-seq']
	meta_df = meta_df[meta_df['Biosample organism'] == 'Mus musculus']
	meta_df = meta_df[meta_df['Experiment target'] == 'CTCF-mouse']
	meta_df = meta_df[meta_df['Project'] == 'ENCODE']
	meta_df = meta_df[meta_df['File analysis title'].str.contains('ENCODE4', na = False)]
	meta_df = meta_df[meta_df['File type'] == 'bed']
	# I used the function metadata_df.nunique() to find the number of unique values in each column, and that helped me figure the necessary columns that contains different values for different dataset that I will download
	meta_df = meta_df[['File accession', 'Output type', 'Experiment accession', 'Assay', 'Donor(s)', 'Biosample term id', 'Biosample term name', 'Biosample type', 'Experiment target', 'File download URL', 'File analysis title']]
	print(meta_df['File accession'].unique())
	print(len(meta_df['File accession'].unique()))
	meta_df.to_csv(output_fn, header = True, index = False, sep = '\t')
	return meta_df

def download_data_one_row(row):
	biosample_term = '-'.join(row['Biosample term name'].split())
	biosample_save_fn = '{id}_{term}.bed.gz'.format(id  = row['File accession'], term = biosample_term)
	save_fn = os.path.join(args.download_folder, biosample_save_fn)
	if not os.path.isfile(save_fn):
		response = requests.get(row['File download URL'])
		file = open(save_fn, 'wb')
		file.write(response.content)
		file.close()
	return 

if __name__ == '__main__':
	meta_df = cleanup_metadata(args.metadata_fn, args.output_fn)
	print('Done cleaning up metadata')
	if args.download:
		meta_df.apply(download_data_one_row, axis = 1) # apply function to download data for each row (each data of the cleaned metadata)
	print('Done downloading data')

