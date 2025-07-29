#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=2_cutadapt_minimap2_tclean-%A_%a.o
#SBATCH --error=2_cutadapt_minimap2_tclean-%A_%a.e


##-------------------------------------------------------------------------

echo Job started on:
date -u

# source config and function
module load Miniconda2/4.3.21
source activate lrp

# load config file provided on command line when submitting job
# Check if a config file was provided on the command line
if [ -z "$1" ]; then
    echo "Error: No config file provided."
    exit 1
fi
config=$(realpath "$1")
echo "Loading config file for project: ${config}" 
source ${config}
source ${SCRIPT_ROOT}/processing/01_source_functions.sh


##-------------------------------------------------------------------------

# create folders
if [ DEMULTIPLEX != TRUE ]; then
  mkdir -p ${WKD_ROOT}/1_basecalled
else
  mkdir -p ${WKD_ROOT}/1_demultiplex ${WKD_ROOT}/1b_demultiplex_merged
fi 
mkdir -p ${WKD_ROOT}/2_cutadapt_merge ${WKD_ROOT}/3_minimap ${WKD_ROOT}/4_tclean ${WKD_ROOT}/5_cupcake
mkdir -p ${WKD_ROOT}/5_cupcake/5_align  
mkdir -p $WKD_ROOT/5_cupcake/5_align/combined
mkdir -p $WKD_ROOT/5_cupcake/6_collapse $WKD_ROOT/5_cupcake/7_sqanti3


##-------------------------------------------------------------------------

if [ "${MULTIPLEXING}" == TRUE ]; then 

   echo "Processing multiple samples one one flow cell"
   
   sample=${ALL_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
   echo ${sample}
   
   # merge each sample into one fastq file 
   if [ "${MERGE_FASTQ}" == TRUE ]; then
       merge_fastq_across_samples ${sample} ${WKD_ROOT}/1_demultiplex ${WKD_ROOT}/1b_demultiplex_merged
   fi

elif [ "${MULTIPLEXING}" != TRUE ]; then

   # merge each sample into one fastq file 
   if [ "${MERGE_FASTQ}" == TRUE ]; then
      for ((n = 0; n <= numSamples -1; n++)); do
        merge_fastq_across_samples ${ALL_SAMPLES_NAMES[$n]} ${RAW_ROOT_DIR[$n]} ${WKD_ROOT}/1_basecalled
      done      
   fi

   if [ "${numSamples}" == "1" ]; then
   
     echo "Processing one sample on one flow cell"
     
     split_merged_fastq ${WKD_ROOT}/1_basecalled
     
     parts=(00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20)
     part=${parts[${SLURM_ARRAY_TASK_ID}]}
     sample=${ALL_SAMPLES_NAMES[0]}_${part}

   else
    
     echo "Processing one sample on one flow cell, but across multiple flow cells"
     
     split_merged_fastq ${WKD_ROOT}/1_basecalled
  
     parts=(00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20)
     
     for sample_name in "${ALL_SAMPLES_NAMES[@]}"; do
       echo ${sample_name}     
       for part in "${parts[@]}"; do       
        samplePart+=("${sample_name}_${part}")     
       done   
     done
     
     sample=${samplePart[${SLURM_ARRAY_TASK_ID}]}
  
    fi 

else

   # Exit with error if none of the conditions are met
   echo "Error: Invalid configuration for demultiplexing, check wiki for combinations."
   exit 1

fi



##-------------------------------------------------------------------------

if [ "${MULTIPLEXING}" == TRUE ]; then 
  
  if [ "${DEMULTIPLEX}" == "FALSE" ] && [ "${DEMULTIPLEX_SOFTWARE}" == "Pychopper" ]; then
      # Run pychopper if DEMULTIPLEX is FALSE and software is Pychopper
      run_pychopper ${WKD_ROOT}/1b_demultiplex_merged/${sample}_merged.fastq ${WKD_ROOT}/2_cutadapt_merge
  
  elif [ "${DEMULTIPLEX_SOFTWARE}" == "Porechop" ] || [ "${SEQUENCING}" == "targeted" ]; then
      # Run post-processing for Porechop or targeted sequencing
      post_porechop_run_cutadapt ${WKD_ROOT}/1b_demultiplex_merged/${sample}_merged.fastq ${WKD_ROOT}/2_cutadapt_merge
  
  elif [ "${ERCC}" == "TRUE" ]; then
      # create a symlink between $WKD_ROOT/1_demultiplex and already demuxed folder (overwrites)
      ln -sfn ${GENOME_WKD_ROOT}/2_cutadapt_merge/* "${WKD_ROOT}/2_cutadapt_merge/"
      
  else
      # Exit with error if none of the conditions are met
      echo "Error: Invalid configuration for demultiplexing, check wiki for combinations."
      exit 1
  fi

else 
  
  run_pychopper ${WKD_ROOT}/1_basecalled/${sample}_merged.fastq ${WKD_ROOT}/2_cutadapt_merge
  
  # create stats output for QC downstream
  PrefixOriginal=${ALL_SAMPLES_NAMES[0]}
  seqkit stats -a ${WKD_ROOT}/1_basecalled/${PrefixOriginal}_merged.fastq > ${WKD_ROOT}/1_basecalled/${NAME}_readstats.txt
fi
  
# map combined fasta to reference genome
run_minimap2 ${WKD_ROOT}/2_cutadapt_merge/${sample}_merged_combined.fastq ${WKD_ROOT}/3_minimap

# run transcript clean on aligned reads
run_transcriptclean ${WKD_ROOT}/3_minimap/${sample}_merged_combined_filtered_sorted.sam ${WKD_ROOT}/4_tclean

# re-align reads
run_pbmm2 ${WKD_ROOT}/4_tclean/${sample}/${sample}_clean.fa ${WKD_ROOT}/5_cupcake/5_align

# filter_alignment <input_name> <input_mapped_dir>
# output = ${sample}_mapped.filtered.bam, ${sample}_mapped.filtered.sorted.bam
filter_alignment ${sample}_mapped ${WKD_ROOT}/5_cupcake/5_align
