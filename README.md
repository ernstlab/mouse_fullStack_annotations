
# full_stack_ChromHMM_annotations for mouse (mm10)
Data of genome annotation from full-stack ChromHMM model trained with 901 datasets assaying 14 chromatin marks in 26 different cell or tissue types of the **mouse** genome. This project is an extension of the published manuscript at <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02572-z"> Genome Biology </a> that introduces an universal (pan-tissue-type) annotation of the **human** genome 
# Download links:
Data of full-stack genome annotations for reference assemblies mm10 can be found <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//mouse_fullStack/annotation_for_publication"> here</a>: 
Within this folder:
- File <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//mouse_fullStack/annotation_for_publication/mm10_100_segments_segments.bed.gz">mm10_100_segments_segments.bed.gz</a> contains a simple four column .bed file of **mouse** full-stack state annotation in **mm10** assembly. The fourth column contains a state label with a prefix number that can be used to order the states. The OverlapEnrichment and NeighborhoodEnrichment commands of ChromHMM with the '-labels' option can compute enrichments for this file and order states based on the prefix number.
- File <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//mouse_fullStack/annotation_for_publication/mm10_100_segments_browser.bed.gz">mm10_100_segments_browser.bed.gz</a> contains a browser file of **mouse** full-stack state annotation in **mm10** assembly. This file is compatible to for UCSC genome browser. Since our training data (901 input data tracks) are in mm10, mm10 is the assembly used for original training and annotation.

- Detailed description of states can be found at <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//mouse_fullStack/state_annotation_processed.csv"> csv file </a>. The excel version of these files, with more results of the states' overlap enrichments with external annotations can be found in our manuscript, which will be made publich soon.

# Folders:
Within each subfolders inside this folder. Currently, our subfolders have not been completely cleaned with documentations. But, we are working on this part of the project. 
- Folder ```relationship_with_ct_spec_state```: contains code to reproduce the plots in Fig. S6 and S7 of the manuscript.


# License:
All code is provided under the MIT Open Acess License
Copyright 2022 Ha Vu and Jason Ernst

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Contact:
We currently trying to comment the code and provide as much details on how to reproduce the results as possible. If you run into problems, please contact Ha Vu (havu73@ucla.edu) 
