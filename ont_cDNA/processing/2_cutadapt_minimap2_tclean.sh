#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
##SBATCH --array=0-52%15 ## uncomment and set number of jobs (number of samples if running script separately)
#SBATCH --output=2_cutadapt_minimap2_tclean-%A_%a.o
#SBATCH --error=2_cutadapt_minimap2_tclean-%A_%a.e


##-------------------------------------------------------------------------

echo Job started on:
date -u

# source config and function
module load Miniconda2/4.3.21
source activate lrp
source ${SCRIPT_ROOT}/processing/01_source_functions.sh

sample=${ALL_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}

##-------------------------------------------------------------------------

# merge each sample into one fastq file 
merge_fastq_across_samples ${sample} ${WKD_ROOT}/1_demultiplex ${WKD_ROOT}/1b_demultiplex_merged

# delinate polyA and polyT sequences, reverse complement polyT sequences, remove polyA from all sequences
post_porechop_run_cutadapt ${WKD_ROOT}/1b_demultiplex_merged/${sample}_merged.fastq ${WKD_ROOT}/2_cutadapt_merge

# map combined fasta to reference genome
run_minimap2 ${WKD_ROOT}/2_cutadapt_merge/${sample}_merged_combined.fasta ${WKD_ROOT}/3_minimap

# run transcript clean on aligned reads
run_transcriptclean ${WKD_ROOT}/3_minimap/${sample}_merged_combined_sorted.sam ${WKD_ROOT}/4_tclean

# re-align reads
pbmm2 align ${WKD_ROOT}/4_tclean/${sample} ${WKD_ROOT}/5_cupcake/5_align

# filter_alignment <input_name> <input_mapped_dir>
# output = ${sample}_mapped.filtered.bam, ${sample}_mapped.filtered.sorted.bam
filter_alignment ${sample}_mapped ${WKD_ROOT}/5_cupcake/5_align
