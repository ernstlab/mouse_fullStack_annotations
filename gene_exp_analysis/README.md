This folder contains code to create figure 2C in the manuscript. Within this folder:

- ```helper.py```: script containing helper functions for file managements. 
- ```get_avg_gene_per_state.py```: this script will calculate the average gene expression per mouse full-stack state, in each available cell type. Please refer to the Methods section "Average gene expression associated with each full-stack state" for full details about how the average expression is calculated. 
- ```avg_gene_exp_per_state_per_ct.txt.gz```: the output of the script ```get_avg_gene_per_state.py```. Inside this file, each row corresponds to a state (raw state names: E1 --> E100), each column corresponds to a sample. This file will then be used to plot figure 2C using the R markdown code ```plot_avg_gene_exp_heatmap.Rmd```
```
usage: get_avg_gene_per_state.py [-h] [--gene_exp_folder GENE_EXP_FOLDER]
                                 [--segment_fn SEGMENT_FN]
                                 [--num_chromHMM_state NUM_CHROMHMM_STATE]
                                 [--output_fn OUTPUT_FN]
                                 [--state_annot_fn STATE_ANNOT_FN]

Calculating avg gene expression in each state, in each available cell
types. The gene_exp_folder should be where we store BingRen's gene
expression data for mouse. segment_fn should be from ChromHMM model for the
mouse

options:
  -h, --help            show this help message and exit
  --gene_exp_folder GENE_EXP_FOLDER
                        Where there gene exp data for different cell types
                        are stored. Each file in this folder is named in
                        the format <ct>-<sample_code>.expr, where the
                        provided gene expression data is in FPKM unit. This
                        data can be downloaded and unzipped from http://chr
                        omosome.sdsc.edu/mouse/download/19-tissues-expr.zip
                        from Shen et al., 2012, 'A map of the cis-
                        regulatory sequences in the mouse genome', Nature.
                        On Ha's system:
                        /data/ENCODE/mouse/gene_exp/19-tissues-expr/.
  --segment_fn SEGMENT_FN
                        segment_fn. NOTE: Because our gene expression data
                        is in mm9, we will need to pass the segment_fn in
                        mm9 here
  --num_chromHMM_state NUM_CHROMHMM_STATE
                        number of chromHMM states
  --output_fn OUTPUT_FN
                        output_fn
  --state_annot_fn STATE_ANNOT_FN
                        state_annot_fn
```
- ```plot_avg_gene_exp_heatmap.Rmd```: Rmarkdown file that will plot figure 2C. 

# License:
All code is provided under the MIT Open Acess License
Copyright 2022 Ha Vu and Jason Ernst

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

