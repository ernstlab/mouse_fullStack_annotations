This folder contains code to obtain enrichment analysis results of mouse full-stack states with CTCF elements profiled in multiple mouse cell types (obtained from ENCODE project). The results are presented in Fig. 2F in the manuscript. 
- First, we obtained the metadata of CTCF experiments from ENCODE project. We downloaded the metadata, recorded in file. ```metadata_from_ENCODE.tsv```.
- Second, we downloaded data of CTCF peaks (bed files) based on download links provided in ```metadata_from_ENCODE.tsv```, using the script ```download_CTCF_data.py```
```
usage: download_CTCF_data.py [-h] --metadata_fn METADATA_FN [--download]
                             --download_folder DOWNLOAD_FOLDER --output_fn
                             OUTPUT_FN

This file will filter out the data needed to get the CTCF peaks in mouse

optional arguments:
  -h, --help            show this help message and exit
  --metadata_fn METADATA_FN
                        metadata file that I got from ENCODE: mouse, CTCF,
                        TF Chip-seq
  --download            If this flag is present, we will download the data
                        of CTCF peaks
  --download_folder DOWNLOAD_FOLDER
                        Where data should be downloaded to
  --output_fn OUTPUT_FN
                        The file of metadata that we will show as
                        additional data for the paper

```
Note: we include a file ```metadata_clean_for_publication.txt``` in this folder, which is actually the output this second step (```download_CTCF_data.py```). But, if you want o replicate completely what we did, you will absolutely need to use the code to downloaded CTCF bed file data from ENCODE.
- Third, you will need to use ChromHMM to do overlap enrichment analysis of the chromatin states with each of the CTCF files, corresponding to CTCF from different cell/tissue types. This is the command that we used to run ChromHMM (version 1.23): 
```
java -jar <full_path_to_chromHMM.jar> -b 1 -lowmem -noimage <full path to segmentation file genome_100_segments.bed.gz> <folder storing CTCF bed files from step 2> <full-path-to output file prefix (without .txt)>
``` 
The output of this step is provided in file ```overlap_ctcf_mm10.txt```.
- Lastly, we calculated the geometric mean and STD of overlap fold enrichment across different cell types. We do not include the annotated code to generate plots here since the procedure is standard (we still include code that works locally but has not been fully clean from local file paths in this folder). If you want to obtain clean, fully commented code for this step, please contact graduate student Ha Vu. 