This folder contains code to do analysis and create plots that replicate Figures S6 and S7 in the manuscript. 
- First, we obtained per-cell-type 15-chromatin state annotations for 66 reference epigenomes/cell types from <a href="https://www.nature.com/articles/s41586-020-2093-3"> Gorkin et al., 2020 </a> , with download links curated and provided by <a href="https://www.nature.com/articles/s41467-021-22653-8#Sec35">Kwon and Ernst, 2021 </a> (in Supplementary File 1, tab g). Note: the provided links show that all bed files are named based on the format ```<ct_name>_15_segments.bed.gz```, and saved into a folder that is given variable name ```ct_segment_folder``` in the provided snakefile ```mouse_random_represent_full_stack.snakefile``` (look into the first few lines of the file if you are confused about what we mean here). Within this folder, we store some example bed files representing per-ct annotations into the folder ```example_perCT_segment```. Also note that the segmentation files from Gorkin et al., 2020 are sorted bed files and the code that we provide here also assumes that the per-ct segmentation files are sorted. If you apply our code to some data other than those provided by Gorkin et al., 2020, you need to sort the bed files files, using the command line ```zcat <bed_file_name> | sort -k1,1 k2,2n > <output_fn> ``` 
- Second, download the file of mouse full-stack chromatin state annotation ```genome_100_segments.bed.gz```, introduced in our manuscript (see the main readme for download links)
- Next, you can start running our snakemake pipeline in ```mouse_random_represent_full_stack.snakefile``` to reproduce the plots in Fig.S6 and Fig. S7. To run the snakemake pipeline to conduct the analysis sequentially (which can take a long time), you can run ```snakemake --cores 1 --snakefile mouse_random_represent_full_stack.snakefile```. To run the snakemake pipeline such that jobs can be run in parallel, you can run ```snakemake -j --cluster "<command to submit jobs on the computing cluster>" --snakefile mouse_random_represent_full_stack.snakefile```, where for our local machine settings, ```<command to submit jobs on the computing cluster>``` is ```qsub -V -l h_rt=4:00:00,h_data=4G```, specifying the time and memory needed for each job in the pipeline to run. 
# Some scripts that are part of the pipeline (we already incorporated these files into the snakemake pipeline. Therefore, we won't go into details about where these scripts are used within the pipeline here)
- File ```helper.py``` contains some useful functions used by most other scripts.
- File ```get_mm_enrichments_across_ct.py```
```
	usage: get_mm_enrichments_across_ct.py [-h] [--input_folder INPUT_FOLDER] [--output_fn OUTPUT_FN] [--num_state_ct_model NUM_STATE_CT_MODEL] [--full_state_annot_fn FULL_STATE_ANNOT_FN]
	                                       [--ct_state_annot_fn CT_STATE_ANNOT_FN]
	Create an excel of the top ct-spec states enriched in each full-stack state
	optional arguments:
	  -h, --help            show this help message and exit
	  --input_folder INPUT_FOLDER
	                        where there are multiple subfolders, each containing enrichment data for different cell type specific model
	  --output_fn OUTPUT_FN
	                        Where the files showing max, mean, min fold enrichment across all ct-state states are reported for all cell types for each state is reported
	  --num_state_ct_model NUM_STATE_CT_MODEL
	                        number of states in the ct-spec models
	  --full_state_annot_fn FULL_STATE_ANNOT_FN
	                        file of the full-stack state annotations
	  --ct_state_annot_fn CT_STATE_ANNOT_FN
	                        file of the ct-spec state annotations
```
- File ```get_ranked_max_enriched_25_state.py```
```
	usage: get_ranked_max_enriched_25_state.py [-h] [--input_folder INPUT_FOLDER] [--output_fn OUTPUT_FN] [--num_state_ct_model NUM_STATE_CT_MODEL] [--full_state_annot_fn FULL_STATE_ANNOT_FN]
	                                           [--ct_state_annot_fn CT_STATE_ANNOT_FN] [--ct_group_fn CT_GROUP_FN]

	Create an excel of the top ct-spec states enriched in each full-stack state

	optional arguments:
	  -h, --help            show this help message and exit
	  --input_folder INPUT_FOLDER
	                        where there are multiple subfolders, each containing enrichment data for different cell type specific model
	  --output_fn OUTPUT_FN
	                        Where the ct-state with maximum fold enrichment across all ct-state states are reported for all cell types for each state is reported
	  --num_state_ct_model NUM_STATE_CT_MODEL
	                        number of states in the ct-spec models
	  --full_state_annot_fn FULL_STATE_ANNOT_FN
	                        file of the full-stack state annotations
	  --ct_state_annot_fn CT_STATE_ANNOT_FN
	                        file of the ct-spec state annotations
	  --ct_group_fn CT_GROUP_FN
	                        file of the annotations of the cell types
```
- File ```calculate_summary_sample_regions.py```
```
	usage: calculate_summary_sample_regions.py [-h] [--all_seed_folder ALL_SEED_FOLDER] [--output_folder OUTPUT_FOLDER] [--num_state_ct_model NUM_STATE_CT_MODEL] [--full_state_annot_fn FULL_STATE_ANNOT_FN]
	                                           [--ct_state_annot_fn CT_STATE_ANNOT_FN] [--ct_group_fn CT_GROUP_FN]

	Create an excel of the top ct-spec states enriched in each full-stack state

	optional arguments:
	  -h, --help            show this help message and exit
	  --all_seed_folder ALL_SEED_FOLDER
	                        where there are multiple subfolders, each containing enrichment data for different cell type specific model
	  --output_folder OUTPUT_FOLDER
	                        Where the ct-state with maximum fold enrichment across all ct-state states are reported for all cell types for each state is reported
	  --num_state_ct_model NUM_STATE_CT_MODEL
	                        number of states in the ct-spec models
	  --full_state_annot_fn FULL_STATE_ANNOT_FN
	                        file of the full-stack state annotations
	  --ct_state_annot_fn CT_STATE_ANNOT_FN
	                        file of the ct-spec state annotations
	  --ct_group_fn CT_GROUP_FN
	                        file of the annotations of the cell types
```