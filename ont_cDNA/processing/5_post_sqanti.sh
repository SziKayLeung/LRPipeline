#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=5_post_sqanti-%A_%a.o
#SBATCH --error=5_post_sqanti-%A_%a.e


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

export dir=$WKD_ROOT/5_cupcake
mkdir -p ${dir}/8_characterise

##-------------------------------------------------------------------------
## Characterisation with CPAT and colour by abundance

# LOGEN: subset cupcake classification file 
# merge cupcake classification file with abundance
# filter cupcake classification file with minimum number of reads and counts
subset_quantify_filter_tgenes.R \
--classfile ${dir}/7_sqanti3/${NAME}_collapsed_RulesFilter_result_classification.txt \
--expression ${dir}/6_collapse/demux_fl_count.csv \
--filter --nsample=${nsamples} --nreads=${nreads} --monoexonic=${monoexonic} --target_genes=${tgenesFile}

# working variables
finalanno=${dir}/7_sqanti3/${NAME}_collapsed_RulesFilter_result_classification.counts_filtered.txt 
finaliso=${dir}/7_sqanti3/${NAME}_collapsed_RulesFilter_result_classification.filtered_isoforms.txt

# LOGEN: subset fasta and gtf using the finalised list of target gene isoforms
subset_fasta_gtf.py --gtf ${dir}/7_sqanti3/${NAME}_collapsed.filtered.gtf -i ${finaliso} -o counts_filtered
subset_fasta_gtf.py --fa ${dir}/7_sqanti3/${NAME}_collapsed_corrected.fasta -i ${finaliso} -o counts_filtered

# run_cpat <input_fasta> <output_name> <output_dir>
run_cpat ${dir}/7_sqanti3/${NAME}_collapsed_corrected_counts_filtered.fa ${NAME} ${dir}/8_characterise

# extract_best_orf <sample> <root_dir>
extract_best_orf ${NAME} ${dir}/8_characterise
