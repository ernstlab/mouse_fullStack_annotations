# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

library(tidyverse)
library(tidyr)
library(dplyr)
library(pheatmap)
library(ggplot2)
biosample_color_dict <- c('adipose'= '#AF5B39', 'bone marrow'= '#375623', 'brain'= '#C5912B', 'embryo'= '#0070C0', 'esc'= '#924965', 'fibroblast'= '#C075C3', 'heart'= '#D56F80', 'intestine'= '#D0A39B', 'kidney'= '#F4B084', 'limb'= '#F182BC', 'liver'= '#9BC2E6', 'lung'= '#E41A1C', 'muscle'= '#F9B6CF', 'neuron'= '#FFD924', 'other'= '#999999', 'placenta'= '#69608A', 'red blood cells'= '#55A354', 'retina'= '#F8CBAD', 'spleen'= '#678C69', 'stomach'= '#E5BDB5', 'testis'= '#ACB9CA', 'thymus'= '#C6E0B4', 'white blood cells'= '#678C69')
assay_color_dict <- c('ATAC'= '#E8F484', 'DNase'= '#DBE680', 'H3K27ac'= '#F7CB4D', 'H3K27me3'= '#A6A6A4', 'H3K36me3'= '#49AC5E', 'H3K4me1'= '#EDF732', 'H3K4me2'= '#F8CBAD', 'H3K4me3'= '#F13712', 'H3K79me1'= '#9DCEA8', 'H3K79me2'= '#377A45', 'H3K9ac'= '#F2A626', 'H3K9me3'= '#677BF6', 'H3ac'= '#E5EAE7', 'H3K79me3'= '#95B77E')
state_group_color_dict <- c('HET'= '#b19cd9', 'quiescent'= '#ffffff', 'enhancers'= '#FFA500', 'WkEnh'= '#ffff00', 'znf'= '#7fffd4', 'Tx'= '#006400', 'artifacts'= '#ffffff', 'DNase'= '#fff44f', 'ATAC'= '#f7f3b5', 'Prom'= '#ff4500', 'ReprPC and DNase'= '#d1cf90', 'DNase'= '#fff44f', 'exon'= '#3cb371', 'BivProm'= '#6a0dad', 'acetylations'= '#fffacd', 'ct_spec enhancers'= '#FFA500', 'ReprPC'= '#C0C0C0', 'weak promoters'= '#800080', 'TxEnh'= '#ADFF2F', 'TSS'= '#FF0000', 'others'= '#fff5ee', 'TxWk'= '#228B22', 'TxEx'= '#9BBB59')
num_state <- 100

read_state_annot_df <- function(state_annot_fn, states_to_plot){
	annot_df <- as.data.frame(read.csv(state_annot_fn, header = TRUE, stringsAsFactors = FALSE, sep = '\t'))
	tryCatch({
		annot_df <- annot_df %>% rename('state_group' = 'group')
		}, error = function (e) {message("tried to change column names in annot_df and nothing worth worrying happened")})
	annot_df <- annot_df[annot_df$state %in% states_to_plot,]
	annot_df <- annot_df %>% arrange(state_order_by_group) # order rows based on the index that we get from state_order_by_group column
	return(annot_df)	
}

draw_emission_subsetted_states <- function(emission_fn, annot_fn, save_fn, states_to_plot){
	emission_df <- as.data.frame(read.csv(emission_fn, sep = '\t', header = TRUE)) %>% arrange(assay_big_group, assay, biosample_big_group, biosample_group, biosample_name)
	colnames(emission_df) <- c('experiments', paste0('S', seq(1, num_state)), colnames(emission_df)[102:length(colnames(emission_df))])
	annot_df <- read_state_annot_df(annot_fn, states_to_plot) # this step will choose only the states that we want to plot
	#### getting the plot_df: columns : experiments, rows: states
	plot_df <- emission_df %>% select(-c('assay', "biosample", "mark", "mark_color", "assay_big_group", "biosample_name", "biosample_group", "biosample_big_group", "biosample_color")) %>% select(c('experiments', paste0('S', annot_df$state))) # choosing only experiments and the states, which are ordered based on the group, all based on the annot_state_df
	plot_df <- as.data.frame(t(plot_df))
	colnames(plot_df) <- plot_df[1,] # experiments
	plot_df <- plot_df[-1,] # get rid of the of experiments row
	rownames(plot_df) <- annot_df$mneumonics
	plot_df <- plot_df %>% mutate_all(as.numeric)
	#######
	####### getting the annot_exp_df for experiments ######
	annot_exp_df <- emission_df %>% select(c('biosample_group', 'assay')) 
	rownames(annot_exp_df) <- emission_df$experiments
	############
	####### getting the annot_state_df for states ####
	annot_state_df <- annot_df %>% select(c('state_group'))
	rownames(annot_state_df) <- annot_df$mneumonics
	############
	p <- pheatmap(plot_df, fontsize = 5, annotation_col = annot_exp_df, annotation_row = annot_state_df, annotation_colors = list(biosample_group = biosample_color_dict, assay = assay_color_dict, state_group = state_group_color_dict), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, fontsize_col = 3, angle_col = 90 , cellheight = 5, filename = save_fn)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3)
{
	stop("wrong command line argument format", call.=FALSE)
}
emission_fn <- args[1] # output from prepare_emission_for_plotting.py
annot_fn <- args[2] # where annotations of states, the state group and state index by groups are stored
#annot_fn <- './data/state_annotation_processed.csv'
#emission_fn <- './data/emissions/emissions_100_for_pheatmap.txt'
#save_fn <- './data/emissions/emissions_100.png'
save_fn <- args[3] # where the figures should be stored
states_to_plot <- as.integer(args[4:length(args)])
print(states_to_plot)
print ('Done get_emission_df_from_chromHMM_emission')
draw_emission_subsetted_states(emission_fn, annot_fn, save_fn, states_to_plot)