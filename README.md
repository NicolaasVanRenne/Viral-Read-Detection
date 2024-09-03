# Viral-Read-Detection
using INVADEseq to detect viral reads from scRNA-seq data


### pathseq_process.pbs: 
This will process subfolders of the raw_data_folder where cellranger output should be located. Output will take the names of these subfolders.
This script was adapted from cell_culture_samples_GEX_pipeline.sh and uses INVADEseq.py.

These files were retrieved from the INVADEseq GitHub: https://github.com/FredHutch/Galeano-Nino-Bullman-Intratumoral-Microbiota_2022/


The pathseq reference database can be found on https://console.cloud.google.com/storage/browser/gatk-best-practices/pathseq/resources

### INVADEseq.py
this code is used when running pathseq_process.pbs

### scRNAseq_pipeline_SeuratV4_Narmada2024_PATHSEQ.R
Code to reproduce the seurat object. The processed seurat object can be downloaded from https://zenodo.org/records/13643041


# Reference
When using this code, please cite the INVADEseq authors: Galeano-Niño et al., Nature 2022
