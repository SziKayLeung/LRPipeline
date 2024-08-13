#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
##SBATCH --array=0-8418%50 ## uncomment and set number of jobs (number of fastq files if running script separately)
#SBATCH --output=1_demux_porechop-%A_%a.o
#SBATCH --error=1_demux_porechop-%A_%a.e


##-------------------------------------------------------------------------

echo Job started on:
date -u

# source function script
module load Miniconda2/4.3.21
source activate lrp
source ${SCRIPT_ROOT}/processing/01_source_functions.sh


##-------------------------------------------------------------------------

SamplePath=${raw_merged_fastq_files[${SLURM_ARRAY_TASK_ID}]}
Sample=$(basename ${SamplePath} .fastq.gz)

echo "Processing ${Sample}"

# 3) run_porechop <raw.fastq.gz> <output_dir>
run_porechop ${SamplePath} ${WKD_ROOT}/1_demultiplex/${Sample} > ${WKD_ROOT}/1b_demultiplex_merged/log/${Sample}.log


##-------------------------------------------------------------------------

# print end date and time
source deactivate
echo Job finished on:
date -u