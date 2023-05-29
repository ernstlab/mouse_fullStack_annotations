# Copyright 2022 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import pandas as pd 
import numpy as np 
import helper 
import os 
import argparse
import sys
sys.path.append('../emission_and_metadata') # full path on Ha's system: /u/home/h/havu73/project-ernst/source/mm10_annotations/emission_and_metadata
import prepare_metadata_files as meta
################## READING INPUT COMMAND LINE ARGS ##################
parser = argparse.ArgumentParser(description = 'This file aims at producing the supplementary file listing the input file metadata and download links.')
parser.add_argument('--ENCODE_meta_fn', type = str, required = False, default = './metadata_controltype_june_2021.tsv', help = 'The file where we can get metadata of experiments from ENCODE')
parser.add_argument('--output_fn', type = str, required = True, help = 'The excel output file. Multiple tabs will be in this file. The output file is part of Additional File 2 in the published paper')
args = parser.parse_args()
helper.check_file_exist(args.ENCODE_meta_fn)
helper.create_folder_for_file(args.output_fn)
#########################################################################################
def read_filter_metadata():
	meta_df = pd.read_csv(args.ENCODE_meta_fn, header = 0, index_col = None, sep = '\t', comment = '#', low_memory = False)
	# meta_df = meta_df[meta_df['File analysis title'].str.contains('ENCODE4', na = False)] # note we need to comment out this line because the read files for different control experiments are actually not part of ENCODE4, and are all blank on this column
	ALL_ASSAY_LIST = ['Control ChIP-seq', 'Histone ChIP-seq', 'DNase-seq', 'ATAC-seq']
	meta_df = meta_df[meta_df['Assay'].isin(ALL_ASSAY_LIST)] # filter only lines whose assay is in ALL_ASSAY_LIST
	COLUMNS_TO_PARSE = ['File accession', 'Output type', 'Experiment accession', 'Assay', 'Biosample term name', 'Experiment target', 'File assembly', 'Run type', 'Derived from', 'Controlled by', 'File analysis title']
	meta_df = meta_df[COLUMNS_TO_PARSE]
	return meta_df

def combine_with_organ_meta_df(file_df, organ_meta_df):
	'''
	file_df can be from functions produce_cell_mark_table_from_chipBamDf or from get_openC_df
	'''
	file_df = file_df.merge(organ_meta_df, how = 'left', left_on = 'Experiment accession', right_on = 'experimentID')
	file_df = file_df[['exp_name', 'Experiment target', 'Biosample term name', 'organ_group', 'cell_group', 'development_group', 'system_group', 'Experiment accession', 'download_link', 'File assembly', 'Run type', 'ctrl_fileCode', 'ctrl_download']]
	return file_df

def check_similar_file_acc_code_with_jason_file(df,  jason_fn):
	jason_df = pd.read_csv(jason_fn, header = None, sep = '\t')
	diff_list = np.setdiff1d(df['File accession'], jason_df[0])
	assert len(diff_list) == 0, 'The expected and observed list of file accession code from jason file {} IS NOT EQUAL.'.format(jason_fn)
	return  

def get_control_df(meta_df):
	# first filter so that we only get the Control Chip-seq experiments
	control_df = meta_df[meta_df['Assay'] == 'Control ChIP-seq']
	bam_control_df = control_df[(control_df['Output type'] == 'alignments') & (control_df['File analysis title'].str.contains('ENCODE4', na = False))] # these are control bam files
	jason_bamControl_fn = os.path.join(ENCODE_meta_fn, 'controlmeta_alignments.txt')
	check_similar_file_acc_code_with_jason_file(bam_control_df,  jason_bamControl_fn)
	# assay is always Control ChIP-seq
	# Biosample term name corresponds to different cell/ tissue types
	# Experiment target is all blank
	# Run type is all blank
	# Controlled by is all blank
	# File analysis title all contains ENCODE4
	read_control_df = control_df[control_df['Output type'] == 'reads'] # these are read files for control samples
	jason_readControl_fn = os.path.join(ENCODE_meta_fn, 'controlmeta_reads.txt')
	check_similar_file_acc_code_with_jason_file(read_control_df, jason_readControl_fn)
	return bam_control_df, read_control_df

def get_chip_df(meta_df):
	bam_chip_df = meta_df[(meta_df['Assay'] == 'Histone ChIP-seq') & (meta_df['Output type'] == 'alignments') & (meta_df['File analysis title'].str.contains('ENCODE4', na = False))] # this will result in 1397 bam files belonging to 704 unique ChIP seq experiments
	jason_bamChip_fn = os.path.join(ENCODE_meta_fn, 'histonemeta_alignments.txt')
	check_similar_file_acc_code_with_jason_file(bam_chip_df, jason_bamChip_fn)
	read_chip_df = meta_df[(meta_df['Assay'] == 'Histone ChIP-seq') & (meta_df['Output type'] == 'reads')]
	jason_readChip_fn = os.path.join(ENCODE_meta_fn, 'histonemeta_reads.txt')
	check_similar_file_acc_code_with_jason_file(read_chip_df,  jason_readChip_fn)	
	return bam_chip_df, read_chip_df

def get_openC_df(meta_df, organ_meta_df):
	'''
	Extract bam files corresponding DNase and ATAC-seq output
	'''
	bam_DNase_df = meta_df[(meta_df['Assay'] == 'DNase-seq') & (meta_df['Output type'] == 'alignments') & (meta_df['File analysis title'].str.contains('ENCODE4', na = False))] # based on code in ENCODE_meta_fn/getmeta.sh
	jason_dnaseBam_fn = os.path.join(ENCODE_meta_fn, 'dnasemeta_alignments.txt')
	check_similar_file_acc_code_with_jason_file(bam_DNase_df, jason_dnaseBam_fn)
	bam_ATAC_df = meta_df[(meta_df['Assay'] == 'ATAC-seq') & (meta_df['Output type'] == 'alignments') & (meta_df['File analysis title'].str.contains('ENCODE4', na = False))] # based on code in ENCODE_meta_fn/getmeta.sh
	jason_ATACBam_fn = os.path.join(ENCODE_meta_fn, 'atacmeta_alignments.txt')
	check_similar_file_acc_code_with_jason_file(bam_ATAC_df, jason_ATACBam_fn)
	openc_df = pd.concat([bam_DNase_df, bam_ATAC_df], ignore_index = True)
	openc_df.loc[:, 'Experiment target'] = openc_df['Assay'].apply(lambda x: x.split('-seq')[0])
	openc_df = openc_df.sort_values(['Experiment target', 'Biosample term name', 'Experiment accession'])
	openc_df['download_link'] = openc_df['File accession'].apply(lambda x: 'https://www.encodeproject.org/files/{c}/@@download/{c}.bam'.format(c = x))
	openc_df['exp_name'] = openc_df['Biosample term name'] + '_' + openc_df['Experiment target'] + '_' + openc_df['Experiment accession']
	openc_df['ctrl_fileCode'] = 'uniform'
	openc_df['ctrl_download'] = ''
	openc_df = openc_df[['exp_name', 'Experiment target', 'Biosample term name', 'Experiment accession', 'download_link', 'File assembly', 'Run type', 'ctrl_fileCode', 'ctrl_download']]
	openc_df.reset_index(inplace = True, drop = True)
	openc_df = combine_with_organ_meta_df(openc_df, organ_meta_df)
	return openc_df


def extract_file_code_one_row(entry, check_set):
	'''
	For each entry of /files/ENCFF001JZL/, /files/ENCFF309GLL/, extrack the file accession code
	'''
	if pd.isnull(entry):
		return []
	code_list = entry.split(', ') # from /files/ENCFF001JZL/, /files/ENCFF309GLL/ to ['/files/ENCFF001JZL/', '/files/ENCFF309GLL/']
	code_list = list(map(lambda x: x[7:-1], code_list)) # to ['ENCFF001JZL', 'ENCFF309GLL']
	if check_set != None:
		results = []
		for code in code_list:
			if code in check_set:
				results.append(code)
		return results
	return code_list

def map_controlRead_to_controlBam(bam_control_df, read_control_df):
	# given the table showing the metadata of control bam and control read files, we will output a dictionary with keys: controlRead --> values list of controlBam that were derived from the control read key
	control_readToBam = {}
	set_conRead = set(read_control_df['File accession']) # a set of files that correspond to control read files
	bam_control_df.loc[:, 'read_file_list'] =  bam_control_df['Derived from'].apply(lambda x: extract_file_code_one_row(x, set_conRead))
	for index, row in bam_control_df.iterrows():
		read_control_list = row['read_file_list']
		for file in read_control_list:
			if not (file in control_readToBam):
				control_readToBam[file] = []
			(control_readToBam[file]).append(row['File accession']) # add the control bam file accession code to the list of files derived from the control read file
	return control_readToBam	

def map_file_accession_to_fileCode_list(first_df, second_df, colname_to_extract):
	'''
	first_df should have columns file accession code, and column Derived from or Controlled by (colname_to_extract). These columns should contain data of the files that the file  in 'File accession' are derivied from or controlled by.
	second_df should have column 'File accession' that contains all the legal files to extract files of columns 'Derived from' and 'Controlled by'
	'''
	set_code_list = set(second_df['File accession']) # set of chip read files
	first_df.loc[:, 'read_file_list'] = first_df[colname_to_extract].apply(lambda x: extract_file_code_one_row(x, set_code_list))
	return dict(zip(first_df['File accession'], first_df['read_file_list']))

def map_bamChip_to_bamCon(chipBam_to_chipRead, chipRead_to_conRead, conRead_to_conBam):
	'''
	mappings from chipBam to a list of controlBam, given the 3 exiting mappings
	'''
	results = {} # keys: chipBam, values: list of conBam
	for chipBam in chipBam_to_chipRead:
		conBam_set = set([])
		chipRead_list = chipBam_to_chipRead[chipBam]
		for chipRead in chipRead_list:
			conRead_list = chipRead_to_conRead[chipRead]
			for conRead in conRead_list:
				conBam_list = conRead_to_conBam[conRead]
				conBam_set.update(conBam_list)
		results[chipBam] = list(conBam_set)
	return results


def match_conBam_to_chipBam(chipBam, chipBam_to_conBam):
	'''
	chipBam is a file accession code
	'''
	try:
		return chipBam_to_conBam[chipBam]
	except:
		return []

def produce_cell_mark_table_from_chipBamDf(bam_chip_df, chipBam_to_conBam, organ_meta_df):
	'''
	chipBam_to_conBam: a dictionary of keys: chipBam, values: list of conBam associated with each chipBam
	'''
	bam_chip_df = bam_chip_df.sort_values(['Experiment target', 'Biosample term name', 'Experiment accession'])
	bam_chip_df.loc[:, 'exp_name'] = bam_chip_df['Biosample term name'] + '_' + bam_chip_df['Experiment target'] + '_' + bam_chip_df['Experiment accession']
	bam_chip_df['control_bam_list'] = bam_chip_df['File accession'].apply(lambda x: match_conBam_to_chipBam(x, chipBam_to_conBam))
	result_df = pd.DataFrame(columns = ['exp_name', 'Experiment target', 'Biosample term name', 'Experiment accession', 'download_link', 'File assembly', 'Run type', 'ctrl_fileCode', 'ctrl_download'])
	for index, row in bam_chip_df.iterrows():
		chip_download = 'https://www.encodeproject.org/files/{c}/@@download/{c}.bam'.format(c = row['File accession'])
		chip_report_row = pd.Series([row['exp_name'], row['Experiment target'], row['Biosample term name'], row['Experiment accession'], chip_download, row['File assembly'], row['Run type']], index = ['exp_name', 'Experiment target', 'Biosample term name', 'Experiment accession', 'download_link', 'File assembly', 'Run type'])
		control_bam_list = row['control_bam_list']
		for ctrl_bam in control_bam_list:
			ctrl_download =  'https://www.encodeproject.org/files/{c}/@@download/{c}.bam'.format(c = ctrl_bam)
			ctrl_report_row =  pd.Series([ctrl_bam, ctrl_download], index = ['ctrl_fileCode', 'ctrl_download'])
			output_row = pd.concat([chip_report_row, ctrl_report_row])
			result_df.loc[result_df.shape[0], :] = output_row
	result_df = combine_with_organ_meta_df(result_df, organ_meta_df)
	return result_df


def save_output(output_fn, chip_output_df, openc_df):
	writer = pd.ExcelWriter(output_fn, engine = 'xlsxwriter')
	chip_output_df.to_excel(writer, sheet_name = 'Histone ChIP-seq')
	openc_df.to_excel(writer, sheet_name = 'ATAC_DNase')
	writer.save()
	return 

if __name__ =='__main__':
	organ_meta_df = meta.get_metadata_from_jason() # output is a df outlinining the cell types of experiment. columns: 'experimentID', 'organ_group', 'cell_group', 'development_group', 'system_group'
	meta_df = read_filter_metadata()
	bam_control_df, read_control_df = get_control_df(meta_df)
	bam_chip_df, read_chip_df = get_chip_df(meta_df)
	chipBam_to_chipRead = map_file_accession_to_fileCode_list(bam_chip_df, read_chip_df, 'Derived from')
	chipRead_to_conRead = map_file_accession_to_fileCode_list(read_chip_df, read_control_df, 'Controlled by')
	conRead_to_conBam = map_controlRead_to_controlBam(bam_control_df, read_control_df)
	chipBam_to_conBam = map_bamChip_to_bamCon(chipBam_to_chipRead, chipRead_to_conRead, conRead_to_conBam)
	chip_output_df = produce_cell_mark_table_from_chipBamDf(bam_chip_df, chipBam_to_conBam, organ_meta_df)
	openc_df = get_openC_df(meta_df, organ_meta_df)
	save_output(args.output_fn, chip_output_df, openc_df)
