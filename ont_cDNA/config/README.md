## Config file

To run the bioinformatics pipeline, a config file of the file paths and parameters are required:

### Arguments to specify the pipeline

-  `export MULTIPEXING=<TRUE><FALSE>` to specify if de-multiplexing is required
    - `TRUE` if more than one sample was sequenced on one flow cell, or de-multiplexing has already been performed

- `export SEQUENCING=<targeted><whole>` to specify which tool to use for de-multiplexing if `MULTIPLEXING=TRUE`
  - `targeted`: a panel genes were enriched for sequencing, run *porechop*
  - `whole`: whole transcriptome profiling, run *pychopper*  

<br>

```
export NAME=AllBDRTargeted  # This will be the output prefix for all downstream files

# if multiplexing was performed
export MULTIPLEXING=TRUE

# sequencing mode: <targeted> <whole>
export SEQUENCING=targeted
```

### Directory paths and other scripts

- `export SCRIPT_ROOT`: path to the LRPipeline GitHub repository scripts
- `export WKD_ROOT`: path to the root directory for output files generated from this pipeline
- `export META_ROOT`: path to the root directory where barcodes and samples csv files are stored
- `export LOGEN_ROOT`: path to the LOGen GitHub repository scripts

<br>

```
## Output root directory filepath (ensure path exists)
export SCRIPT_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LRPipeline/ont_cDNA
export WKD_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/D_ONT
export META_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/0_metadata/B_ONT

## Generate folder structure in root output directory 
cd ${WKD_ROOT}
mkdir -p 1_demultiplex 1b_demultiplex_merged 2_cutadapt_merge 3_minimap 4_tclean 5_cupcake
mkdir -p ${WKD_ROOT}/1_demultiplex/Batch2 ${WKD_ROOT}/5_cupcake/5_align
mkdir -p $WKD_ROOT/5_cupcake/5_align/combined $WKD_ROOT/5_cupcake/6_collapse $WKD_ROOT/5_cupcake/5_align/combined_fasta $WKD_ROOT/5_cupcake/7_sqanti3

## SKLeung scripts
export LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous # no need to update
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing # no need to update
SUBSETPOLYTAILS=$LOGEN_ROOT/assist_ont_processing/subset_polyA_polyT.py # no need to update
```

### References

- `export GENOME_FASTA`: path to the reference fasta for alignment
- `export GENOME_GTF`: path to the reference gtf for alignment  

<br>

```
## Reference data filepaths
export GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
export GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/annotation/gencode.v40.annotation.gtf
```

### Raw data

- `export numSamples`: number of samples sequenced across all flow cells 
- `export raw_merged_fastq_files`: path to the root directory containing raw passed fastq files
- `export BARCODE_CONFIG`: path to barcode csv
- `export SAMPLE_ID`: path to sample csv

<br>

```
# number of total samples (across all flow cells)
export numSamples=53

## ONT raw data
# sequentially specify paths of raw fastq files from multiple flow cells
export RAW_ROOT_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/sorted_nuclei/RNA/human
export raw_merged_fastq_files=${RAW_ROOT_DIR}/P0075_20230216_10916/adult_sorted_nuclei/20230216_1558_1G_PAK65840_66586f9d/fastq_pass

# export barcode (if multiplexing=TRUE) <barcode number> <sample name>
BARCODE_CONFIG=${SCRIPT_ROOT}/config/barcode.csv
export ALL_SAMPLES_NAMES=($(awk -F "\"*,\"*" '{print $1}' ${BARCODE_CONFIG})) # no need to update

# sample names to replace barcode names downstream 
SAMPLE_ID=${SCRIPT_ROOT}/config/sample_id.csv
```

### Software 

Path to other third-party tools 

<br>

```
# Path to software root directory
export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software

# QC
export MINIONQC=${SOFTDIR}/minion_qc/MinIONQC.R

# Porechop
export PORECHOP=${SOFTDIR}/Porechop/porechop-runner.py

# TranscriptClean
export TCLEAN=${SOFTDIR}/TranscriptClean/TranscriptClean.py

# IsoSeq3 cupcake 
export CUPCAKE=${SOFTDIR}/cDNA_Cupcake
export ANNOTATION=$CUPCAKE/annotation
export SEQUENCE=$CUPCAKE/sequence
export PYTHONPATH=$PYTHONPATH:$SEQUENCE

# SQANTI
export SQANTI3_DIR=${SOFTDIR}/SQANTI3
export SQANTI_JSON=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json
CAGE_PEAK=$SQANTI3_DIR/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
POLYA=$SQANTI3_DIR/data/polyA_motifs/mouse_and_human.polyA_motif.txt
```
