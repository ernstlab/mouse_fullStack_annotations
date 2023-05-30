This folder contains code that produces excel files associated with different overlap enrichment analysis presented in the paper. In these analyses, for each state in the mouse full-stack annotation, we calculate the overlap fold-enrichment between the state and an annotation of interests (for example, with different chromosomes, different genome contexts from RefSeq, different classes of repeat elements). Please refer to the Method subsection "External annotation sources" for a list of external genome annotations that we use to overlap with the mouse full-stack states.
The process of doing overlap enrichment includes (1) Run ```ChromHMM OverlapEnrichment``` to calculate the fold-overlap, and (2) Run ```create_excel_painted_overlap_enrichment.py``` to create excel file that visualize the results.
(1) We run ```ChromHMM OverlapEnrichment``` with <a href="http://compbio.mit.edu/ChromHMM/">ChromHMM v.1.23</a> following the manual's instruction. In particular, this is the format of the command to run ChromHMM:
```
java -jar <path/to/ChromHMM.jar> OverlapEnrichment -b 1 -lowmem -noimage <path/to/genome_100_segments.bed.gz: path to bed file of chromatin state segmentation> <path/to/folder/or/bed/file/of/external/annotations> <path/output/file/and/output/prefix>
```
The ```-b 1``` flag tells ChromHMM to calculate overlap enrichment at 1bp resolution
The ```-lowmem``` and ```-noimage``` flags tell ChromHMM to calculate overlap enrichment using lower memory mode, and without generating output heatmap file.
(2) We visualize the results of overlap enrichment using the following scripts:
- ```helper.py```: script containing helper functions for file managements. 
- ```create_excel_painted_overlap_enrichment.py```: File to create excel files showing the overlap enrichment results
```
python create_excel_painted_overlap_enrichment_output.py <input_fn> <output_fn> <context_prefix> <state_annot_fn>
input_fn: the output of chromHMM OverlapEnrichment showing the fold enrichment between each full-stack state and each external genome annotation
output_fn: the name of the excel file that you want to create
context_prefix: the prefix to all the enrichment contexts in this input_fn. For example, if we do enrichments with gnomad variants of varying maf, the column names in input_fn are 'maf_0_0.1<etc>', we can specify context_prefix to gnomad. If we dont need such prefix, specify to empty string. This is useful when we try to write comment on states that are most enriched in each genomic context.
state_annot_fn: where we get all the data of characterization of states, recommended ../state_annotation_processed.csv
```

# License:
All code is provided under the MIT Open Acess License
Copyright 2022 Ha Vu and Jason Ernst

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

