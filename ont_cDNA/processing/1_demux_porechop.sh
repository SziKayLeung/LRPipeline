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
#SBATCH --array=16
#SBATCH --output=1_demux_porechop-%A_%a.o
#SBATCH --error=1_demux_porechop-%A_%a.e


##-------------------------------------------------------------------------

echo Job started on:
date -u

# source function script
module load Miniconda2/4.3.21
source activate sqanti2_py3
source ${SCRIPT_ROOT}/processing/01_source_functions.sh

# load config file provided on command line when submitting job
# Check if a config file was provided on the command line
if [ -z "$1" ]; then
    echo "Error: No config file provided."
    exit 1
fi
config=$(realpath "$1")
echo "Loading config file for project: ${config}"
source ${config}

##-------------------------------------------------------------------------

splitFastq=($(ls ${WKD_ROOT}/1_demultiplex/split/splitfastq*.txt))
SamplePath=${splitFastq[${SLURM_ARRAY_TASK_ID}]}
Sample=$(basename ${SamplePath} .txt)

echo "Merging fastq from ${SamplePath}"
if [ ! -f ${WKD_ROOT}/1_demultiplex/split/${Sample}.fastq ]; then 
  cat $(grep -v '^#' ${SamplePath}) > ${WKD_ROOT}/1_demultiplex/split/${Sample}.fastq.gz
  gunzip ${WKD_ROOT}/1_demultiplex/split/${Sample}.fastq.gz
fi

echo "Processing ${Sample}"

# 3) run_porechop <raw.fastq.gz> <output_dir>
if [ ! -f ${WKD_ROOT}/1_demultiplex/split/${Sample}.log ]; then
  splitFastqPath=${WKD_ROOT}/1_demultiplex/split/${Sample}.fastq
  run_porechop ${splitFastqPath} ${WKD_ROOT}/1_demultiplex/split/${Sample} > ${WKD_ROOT}/1_demultiplex/split/${Sample}.log
fi


##-------------------------------------------------------------------------

# print end date and time
source deactivate
echo Job finished on:
date -u
