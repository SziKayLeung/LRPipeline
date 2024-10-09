#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --output=3_merged_collapse_sqanti3.o
#SBATCH --error=3_merged_collapse_sqanti3.e


##-------------------------------------------------------------------------

echo Job started on:
date -u

# source config and function
module load Miniconda2/4.3.21
source activate nanopore
source ${SCRIPT_ROOT}/processing/01_source_functions.sh
export dir=$WKD_ROOT/5_cupcake

# load config file provided on command line when submitting job
if [ -z "$1" ]; then
    echo "Error: No config file provided."
    exit 1
fi
echo "Loading config file for project: $1" 
source $1


##-------------------------------------------------------------------------

# replace sample barcodes with sample names
replace_filenames_with_csv.py --copy --ext=filtered.bam -i=$WKD_ROOT/5_cupcake/5_align -f=${SAMPLE_ID} -d=${dir}/5_align/combined 
replace_filenames_with_csv.py --copy --ext=fa -i=$WKD_ROOT/5_cupcake/5_align -f=${SAMPLE_ID}  -d=${dir}/5_align/combined_fasta 

# merge all aligned files
allfilteredmapped=($(ls ${dir}/5_align/combined/*filtered.bam)) 
ls ${allfilteredmapped[@]}
samtools merge -f ${dir}/6_collapse/${NAME}_mapped.filtered.sorted.bam ${allfilteredmapped[@]}

# collapse isoforms
run_isoseq_collapse ${dir}/6_collapse/${NAME}_mapped.filtered.sorted.bam ${NAME} 

# demuliplex_collapsed_isoforms <input_directory_fasta> <input_collapsed_directory> <output_name>
demuliplex_collapsed_isoforms ${dir}/5_align/combined_fasta ${dir}/6_collapse ${NAME}

# sqanti3
run_sqanti3 ${dir}/6_collapse/${NAME}_collapsed.gff ${dir}/7_sqanti3

##-------------------------------------------------------------------------

# print end date and time
source deactivate
echo Job finished on:
date -u