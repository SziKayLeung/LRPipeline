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
#SBATCH --output=run_pipeline.o
#SBATCH --error=run_pipeline.e


##-------------------------------------------------------------------------

echo Job started on:
date -u

# load config file provided on command line when submitting job
config=$(realpath "$1")
echo "Loading config file for project: ${config}" 
source ${config}

if [ "${DEMULTIPLEX}" == "TRUE" ]; then 
  if [ "${SEQUENCING}" == "targeted" ]; then
    
    ls ${raw_merged_fastq_files}/*fastq* > ${WKD_ROOT}/1_demultiplex/split/all_fastq.txt
    split -n l/20 -d --additional-suffix=.txt ${WKD_ROOT}/1_demultiplex/split/all_fastq.txt ${WKD_ROOT}/1_demultiplex/split/splitfastq_
    echo "Performed targeted sequencing or use of custom barcodes: using Porechop for demultiplexing primers and barcodes"
    jobid1=$(sbatch ${SCRIPT_ROOT}/processing/1_demux_porechop.sh ${config} | awk '{print $NF}')
  else
    echo "Performed whole transcriptome sequencing or use of standard barcodes: using Pychopper for demultiplexing primers and barcodes"
    jobid1=$(sbatch ${SCRIPT_ROOT}/processing/1_demux_pychopper.sh)
  fi
else
  # create a symlink between $WKD_ROOT/1_demultiplex and already demuxed folder (overwrites)
  ln -sf ${DEMULTIPLEX_DIR}/* "${WKD_ROOT}/1_demultiplex"
  echo "Demultiplexing already performed"
fi

jobid2=$(sbatch --array=0-$((numSamples - 1)) --dependency=afterok:$jobid1 ${SCRIPT_ROOT}/processing/2_cutadapt_minimap2_tclean.sh ${config} | awk '{print $NF}')

# isoseq-collapse, sqanti3
sbatch --dependency=afterok:$jobid2 ${SCRIPT_ROOT}/processing/3_merged_collapse_sqanti3.sh

# QC
sbatch --dependency=afterok:$jobid2 ${SCRIPT_ROOT}/processing/4_QC.sh ${config}