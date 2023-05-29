## NOTES

This folder provides code that were used to generate Additional File 2, tabs 'trainData_HistoneChip' and 'trainData_ATACDNase' outlining the metadata and download links to the datasets that we used to train the mouse full-stack model. Within this folder:
- ```metada_conntroltype_june_2021.tsv```: tsv file downloaded from ENCODE portal for the Chip-seq and ATAC-seq experiments, including control experiments and their files, that we are using for this study.
- ```extract_cell_mark_table.py```: code to produce Additional File 2, tabs 'trainData_HistoneChip' and 'trainData_ATACDNase' in the published paper. This code will extract the metadata for the experiments that we will use for this paper, and their corresponding control files. The control experiments are needed for binarizing input data so that we can use as input to the ChromHMM LearnModel function (see Methods section of the paper).
```
usage: extract_cell_mark_table.py [-h] [--ENCODE_meta_fn ENCODE_META_FN]
                                  --output_fn OUTPUT_FN

This file aims at producing the supplementary file listing the input file
metadata and download links.

optional arguments:
  -h, --help            show this help message and exit
  --ENCODE_meta_fn ENCODE_META_FN
                        The file where we can get metadata of experiments
                        from ENCODE
  --output_fn OUTPUT_FN
                        The excel output file. Multiple tabs will be in
                        this file. The output file is part of Additional
                        File 2 in the published paper
```
- ```helper.py```: script containing helper functions for file managements. 

## LICENSE
Copyright 2022 Ha Vu (havu73@ucla.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.