#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=7:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=200G # specify bytes memory to reserve


##-------------------------------------------------------------------------

echo Job started on:
date -u

# source config and function
module load Miniconda2/4.3.21

# load config file provided on command line when submitting job
if [ -z "$1" ]; then
    echo "Error: No config file provided."
    exit 1
fi

echo "Loading config file for project: $1" 
source $1

source activate lrp

#rootDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/sorted_nuclei/RNA/human/combined/v2
#LRPipelineDir=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LRPipeline
#LOGenDir=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
#manifestQC=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/sorted_nuclei/RNA/human/combined/v2/0_metadata/manifestQC.csv
#ERCC_WKD_ROOT=NULL
#ERCC_WKD_ROOT="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/sorted_nuclei/RNA/human/combined/v2"
#studyName="Humansortednucleidataset"

if [ $MULTIPLEXING == FALSE ]; then
  mkdir -p $WKD_ROOT/3_minimap/partStats
  mv $WKD_ROOT/3_minimap/*mappedstats* $WKD_ROOT/3_minimap/partStats/
  stats=$(ls $WKD_ROOT/3_minimap/partStats/*mappedstats*)
  cat ${stats[@]} > $WKD_ROOT/3_minimap/${NAME}_merged_combined_mappedstats.txt
  
  mkdir -p $WKD_ROOT/5_cupcake/5_align/partStats
  mv $WKD_ROOT/5_cupcake/5_align/PAF/*mappedstats* $WKD_ROOT/5_cupcake/5_align/partStats/
  stats=$(ls $WKD_ROOT/5_cupcake/5_align/partStats/*mappedstats*)
  cat ${stats[@]} > $WKD_ROOT/5_cupcake/5_align/PAF/${NAME}_merged_combined_mappedstats.txt
  
fi

if [ -f ${WKD_ROOT}/QC_input.RData ]; then 
	echo "Input for QC ready"
else 
	echo "Prepare input for QC"
	Rscript ${SCRIPT_ROOT}/QC/read_QC_files.R -m ${manifest} -r ${WKD_ROOT}
fi 

if [ "${ERCC_WKD_ROOT}" == "NULL" ]; then
  ERCCDirParam="NULL"
else
  ERCCDirParam="'${ERCCDir}'"  # Quote the ERCCDir path
fi

Rscript -e "rmarkdown::render('${SCRIPT_ROOT}/QC/QC_report.Rmd', output_file='${WKD_ROOT}/QC_report.html', 
  params = list(
    rootDir = '${WKD_ROOT}', 
    LRPipelineDir = '${SCRIPT_ROOT}', 
    LOGenDir = '${LOGEN_ROOT}', 
    ERCCDir = ${ERCCDirParam}
  ))"
