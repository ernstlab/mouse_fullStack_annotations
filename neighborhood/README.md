This folder contains code to obtain neighborhood enrichment analysis results of mouse full-stack states, as presented in Fig. 2D-E of the manuscript.
- Files ```genome_100_RefSeqTES_neighborhood.txt``` and ```genome_100_RefSeqTSS_neighborhood.txt``` show the output results of running ChromHMM neighborhood enrichments for the 100 mouse full-stack states with the neighborhood of annotated TSS and TES from RefSeq. The bed files of annotated TSS and TES are both provided along with ChromHMM software. 
- You will need to install some R packages if you have not had them in order to use the next two code files. In your R terminal, use the following commands:
```
install.packages("tidyverse")
install.packages("tidyr")
install.packages("dplyr")
install.packages('ggplot2')
install.packages(pheatmap)
install.packages(reshape2)
```
- File ```draw_2DLine_neighborhood_enrichment.R``` contain code to plot figure 2D-E. You can use this command:
```
Rscript draw_2DLine_neighborhood_enrichment.R genome_100_RefSeqTSS_neighborhood.txt genome_100_RefSeqTSS_neighborhood_linePlot.png genome_100_RefSeqTES_neighborhood.txt genome_100_RefSeqTES_neighborhood_linePlot.png ../state_annotation_processed.csv
```

- File ```draw_neighborhood_enrichment.R``` contains code to plot figures S2 in the manuscript. You can use these two commands:
```
Rscript draw_neighborhood_enrichment.R genome_100_RefSeqTSS_neighborhood.txt genome_100_RefSeqTSS_neighborhood_heatMap.png ../state_annotation_processed.csv

Rscript draw_neighborhood_enrichment.R genome_100_RefSeqTES_neighborhood.txt genome_100_RefSeqTES_neighborhood_heatMap.png ../state_annotation_processed.csv
```
Sorry for the short instructions, it's quite self-explanatory, and let us know if you need any questions.
