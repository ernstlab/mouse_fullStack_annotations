This folder contains code to produce Figure 1 of the manuscript. 

- ```helper.py```: script containing helper functions for file managements. 
- Folder ```data```: contains emission probability file and different metadata files that specify the color codes for different chromatin marks, cell types, etc.
- ```prepare_emission_for_plotting.py```: Code to transform data from ```./data/emissions_100.txt``` to ```./data/emissions_100_for_pheatmap.txt```, with added columns showing the associated chromatin mark and cell types fro each experiment. 
- ```rank_chromMark_cellType_emission.py```: Code to produce excel files that correspond to Fig. 1B-C.
```
python rank_chrom_mark_celltype_emission.py <emission_fn> <output_fn> <num_top_marks>
emission_fn: should be ./data/emissions_100.txt
output_fn: execl file where the output data of ranked cell type and chrom marks are stored
num_top_marks: number of top marks that we want to report, recommended value: 100
```
- ```draw_heatmap_emission_mm10_full_stack.Rmd```: Code to produce Fig. 1A, and individual plots corresponding to indidivual states (not presented in the paper, but these plots are helpful in understanding the states' biological implications).
- ```draw_emission_subsetted_states.R```: Code to draw emission probabilities for an user-specified subset of states. This function did not involve in creating any figures in the paper. It may be helpful to users, but it requires that you call states by their raw numbers (column ```state``` in ```../state_annotation_processed.csv```) 
```
Rscript draw_emission_subsetted_states.R <emission_fn: typically ./data/emissions_100_for_heatmap.txt, outputted from prepare_emission_for_plotting.py> <annot_fn: state annotations, ../state_annotation_processed.csv> <save_fn: where we save the output figure, make sure that the path to this file is already created (folders are created)> <states_to_plot: a space-separated sequence of state numbers (raw numbers, 1--> 100), to plot>
```
# License:
All code is provided under the MIT Open Acess License
Copyright 2022 Ha Vu and Jason Ernst

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

