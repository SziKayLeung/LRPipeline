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
#SBATCH -o /dev/null  
#SBATCH -e /dev/null   

##-------------------------------------------------------------------------

echo Job started on:
date -u

# source config and function
module load Miniconda2/4.3.21

# load config file provided on command line when submitting job
config=$(realpath "$1")
if [ -z "$1" ]; then
    echo "Error: No config file provided."
    exit 1
else
    echo "Loading config file for project: ${config}" 
    source ${config}
    source ${SCRIPT_ROOT}/processing/01_source_functions.sh
fi

##-------------------------------------------------------------------------

# create folders
if [ ${DEMULTIPLEX:-TRUE} == "TRUE" ]; then
  mkdir -p ${WKD_ROOT}/1_basecalled
else
  mkdir -p ${WKD_ROOT}/1_demultiplex ${WKD_ROOT}/1b_demultiplex_merged
fi 
mkdir -p ${WKD_ROOT}/0_log
mkdir -p ${WKD_ROOT}/2_cutadapt_merge ${WKD_ROOT}/3_minimap ${WKD_ROOT}/4_tclean ${WKD_ROOT}/5_cupcake
mkdir -p ${WKD_ROOT}/5_cupcake/5_align  
mkdir -p $WKD_ROOT/5_cupcake/5_align/combined
mkdir -p $WKD_ROOT/5_cupcake/6_collapse $WKD_ROOT/5_cupcake/7_sqanti3

# Redirect output manually to ensure it goes there
exec > >(tee -a "${WKD_ROOT}/0_log/run_pipeline.o") 2> >(tee -a "${WKD_ROOT}/0_log/run_pipeline.e" >&2)


##-------------------------------------------------------------------------

# only if multiplexing is performed 
if [ "${MULTIPLEXING}" == TRUE ]; then 

  if [ "${DEMULTIPLEX}" == TRUE ]; then 

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

fi


if [ "${MULTIPLEXING}" == TRUE ]; then 
  jobid1=$(sbatch --array=0-$((numSamples - 1)) ${SCRIPT_ROOT}/processing/2_cutadapt_minimap2_tclean.sh ${config} | awk '{print $NF}')
else
  totalSplitSamples=$((numSamples * 20))
  jobid1=$(sbatch --array=0-$((totalSplitSamples - 1))%5 ${SCRIPT_ROOT}/processing/2_cutadapt_minimap2_tclean.sh ${config} | awk '{print $NF}')
fi

# isoseq-collapse, sqanti3
sbatch --dependency=afterok:$jobid1 ${SCRIPT_ROOT}/processing/3_merged_collapse_sqanti3.sh ${config}

# QC
#sbatch ${SCRIPT_ROOT}/processing/4_QC.sh ${config}
