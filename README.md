
# full_stack_ChromHMM_annotations for mouse (mm10)
Data of genome annotation from full-stack ChromHMM model trained with 901 datasets assaying 14 chromatin marks in 26 different cell or tissue types of the **mouse** genome. The paper has been published on <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02994-x"> Genome Biology </a>. This project is an extension of the published manuscript at <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02572-z"> Genome Biology </a> that introduces an universal (pan-tissue-type) annotation of the **human** genome
# Download links:
Data of full-stack genome annotations for reference assemblies mm10 can be found <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS/full_stack/full_stack_annotation_public_release/mm10/"> here</a>: 
Within this folder:
- File <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//full_stack/full_stack_annotation_public_release/mm10/mm10_100_segments_segments.bed.gz">mm10_100_segments_segments.bed.gz</a> contains a simple four column .bed file of **mouse** full-stack state annotation in **mm10** assembly. The fourth column contains a state label with a prefix number that can be used to order the states. The OverlapEnrichment and NeighborhoodEnrichment commands of ChromHMM with the '-labels' option can compute enrichments for this file and order states based on the prefix number.
- File <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//full_stack/full_stack_annotation_public_release/mm10/mm10_100_segments_browser.bed.gz">mm10_100_segments_browser.bed.gz</a> contains a browser file of **mouse** full-stack state annotation in **mm10** assembly. This file is compatible to for UCSC genome browser. Since our training data (901 input data tracks) are in mm10, mm10 is the assembly used for original training and annotation.

- Detailed description of states can be found at <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//full_stack/full_stack_annotation_public_release/mm10/state_annotation_processed.tsv"> tsv file </a>. The excel version of these files, with more results of the states' overlap enrichments with external annotations can be found in our paper, Additional Files 3-5.

- State group meanings:
```
mGapArtf - assembly gaps and artifacts
mQuies - quiescent
mHet - heterochromatin
mZNF - Zinc finger genes
mReprPC - polycomb repressed
mReprPC_openC - polycomb repressed and open chromatin
mOpenC - open chromatin
mEnhA - active enhancers
mEnhWk - weak enhancers
mTxEnh - transcribed enhancers
mTx - transcription
mTxEx - transcription and exons
mTxWk - weak transcription
mBivProm - bivalent promoters
mPromF - promoter flank
mTSS - transcription start sites
```

# Track hubs on UCSC genome browser:
You can view the full-stack annotations for the mouse (presented in this manuscript) OR the human (presented in the <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02572-z"> "sister" manuscript</a>), by using the <a href="ttps://public.hoffman2.idre.ucla.edu/ernst/2K9RS//full_stack/full_stack_annotation_public_release/hub.txt"> track hub link</a>. We provide a very detailed step-by-step instruction on how to view the full-stack annotations using the provided track hub link in the tutorial file <a href="https://github.com/ernstlab/mouse_fullStack_annotations/blob/main/view_ucsc_genome_browser.pptx">view_ucsc_genome_browser.pptx</a>.

# Folders:
Within each subfolders inside this folder. Each subfolder contains its own readme file.  
- Folder ```relationship_with_ct_spec_state```: contains code to reproduce the plots in Fig. S6 and S7 of the manuscript.
- Folder ```process_CTCF_data```: contains code to reproduce Fig. 2F of the manuscript.
- Folder ```neighborhood```: contains code to reproduce Fig. 2D-E and S2 of the manuscript.
- Folder ```emission_and_metadata```: contains code to reproduce Fig. 1 of the manuscript
- Folder ```gene_exp_analysis```: contains code to reproduce Fig. 2C of the manuscript
- Folder ```overlap```: contains scripts of calculate the enrichments of mouse full-stack states with different genome annotation. 
- Folder ```learn_model```: contains scripts to reproduce tabs 'trainData_HistoneChip' and 'trainData_ATACDNase' in Additional File 2 of the paper

# License:
All code is provided under the MIT Open Acess License
Copyright 2022 Ha Vu and Jason Ernst

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Contact:
If you run into problems, please contact Ha Vu (havu73@ucla.edu) 
