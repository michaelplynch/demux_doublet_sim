# demux_doublet_sim

Repository for Nextflow pipeline used in demuxSNP demultipelxing paper

## Overall workflow

1. Simulate doublets  
Scripts/templates leveraged from Weber et al. (DOI: https://doi.org/10.1093/gigascience/giab062)
- Start with demultiplexed/individual bams
- Barcode suffix is replace with sample key e.g. K1, K2 etc.
- bam files are merged.
- Lookup file generated with barcodes randomly selected to satisfy required proportion of doublets.
- Barcodes in bam renamed according to lookup to simulate doublets.  
2. Benchmark methods  
Tests demuxSNP (hybrid), HTOreader (hybrid) and souporcell.
- Experiments 1: Vary doublet rate
- Experiment 2: Vary SNP subsetting

## Inputs

Most inputs are specified in nextflow.config:
    container__souporcell: path to souporcell apptainer image, ideally at top level of project.  
    bam_path: Path to demultiplexed bam files.  
    barcodes_path: Path to demultiplexed barcodes.  
    tenx: Path to barcodes.tsv, features.tsv and matrix.mtx files from multiplexed 10X output.  
    common_variants: common variants e.g. from 1K genome project.  
    ref: path to reference genome, ideally in data/input directory.  

Doublet simulation parameters are specified in params_ccrcc.csv
The workflow caters for subsampling (also specified in params_ccrcc.csv) although this was not explored in the paper.

## Outputs

Folder for each simulated scenario (e.g. seed, % doublets, number of genes used to subset)
SingleCellExperiment object in each demuxSNP folder.

## Known issues

Input files used by souporcell/apptainer need to be stored below the image.
Apptainer must be bound to the project directory (variable in nextflow.config).