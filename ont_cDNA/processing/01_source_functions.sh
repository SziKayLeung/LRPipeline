source activate lrp 
export PATH=$PATH:${LOGEN_ROOT}/transcriptome_stats
export PATH=$PATH:${LOGEN_ROOT}/compare_datasets
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
export PATH=$PATH:${LOGEN_ROOT}/merge_characterise_dataset
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
SUBSETPOLYTAILS=$LOGEN_ROOT/assist_ont_processing/subset_polyA_polyT.py

# Get the absolute directory where 01_source_functions.sh is located
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Check if the script directory is inside a Git repository
if git -C "$script_dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    # Get the latest commit hash from the Git repository
    current_commit_hash=$(git -C "$script_dir" rev-parse HEAD)
    echo "LRPipeline latest git commit hash: $current_commit_hash"
else
    echo "Not in a Git repository. Skipping git commit hash check."
    current_commit_hash="unknown"
fi

# 1) run_merge <raw_directory> <sample_output_name>
# output: <sample_output_name>.merged.fastq 
run_merge(){

  echo "Merging following fastq files"
  FASTQ=$(ls $1/*) 
  echo ${FASTQ}

  cat ${FASTQ} > $2
  echo "Merge of Samples successful: output to $2"
  
}

# convertfasta2fastq <input_fastq> <output_fasta>
convertfasta2fastq(){

  if [ ! -f $2 ]; then 
    seqtk seq -A $1 > $2
  fi 
}

merge_fastq_across_samples(){
  # variables 
  gval=$1
  input_dir=$2
  output_dir=$3
    
  if [ -f ${output_dir}/${gval}_merged.fastq ]; then
    echo ${output_dir}/${gval}_merged.fastq
    echo -e "Merging ${gval}: \e[32mCompleted\e[0m"
  
  else
    echo "Merging ${gval}"
    
    if [ "${DEMULTIPLEX_SOFTWARE}" == "Porechop" ]; then
        fastq=$(find "$input_dir" -type f -name "*${gval}*")
    else
        # special characters (like trailing spaces, newline characters, or non-printing characters) needs to be removed to correctly access paths
        clean_path=$(echo "${input_dir}/${gval}" | tr -cd '\11\12\15\40-\176')
        echo "$clean_path"
        fastq=$(ls ${clean_path}/*fa* 2>/dev/null)
    fi
    
    num_files=$(echo "$fastq" | wc -w)
    if [ $num_files -eq 0 ]; then
        echo "Merging failed"
        exit 1
    fi 
    
    echo "Number of files to concatenate: $num_files"
    echo "$fastq" > ${output_dir}/${gval}_file_list.txt
    
    # Check if the files are gzipped or plain fastq
    if echo "$fastq" | grep -q ".gz$"; then
      # If files are gzipped, concatenate and output as gzipped
      echo "Concatenating gzipped files and unzip..."
      zcat $fastq > ${output_dir}/${gval}_merged.fastq
    else
      # If files are not gzipped, concatenate as plain fastq
      echo "Concatenating fastq files..."
      cat $fastq > ${output_dir}/${gval}_merged.fastq
    fi
    
    seqkit stats -a ${output_dir}/${gval}_merged.fastq > ${output_dir}/${gval}_readstats.txt
  fi 
  
}

# 2) run_QC <sample> <sequencing_summary> <bam_input> <output_dir>
run_QC(){
    # variables
    sample=$1
    sequencing_summary=$2
    bam_input=$3
    output_dir=$4

    echo "Processing: $1"
    cd $output_dir
    pycoQC --summary_file $sequencing_summary --bam_file $bam_input -o $sample"_QC.html"
    Rscript ${MINIONQC} -i $sequencing_summary -s TRUE -o $output_dir

}

# 3) run_pychopper <input_fastq> <output_dir>
# input: <demultiplexed_merged>.fastq 
# output: <output_directory>/<sample>_merged_combined.fasta
run_pychopper(){

  name=$(basename $1 .fastq)
  
  echo "Processing ${name} for pychopper" 
  
  if [ -f $2/${name}_combined.fasta ]; then
  
    echo -e "Pychopper: \e[32mCompleted\e[0m"
  
  else

    pychopper -r $2/${name}_pychopperReport.pdf $1 $2/${sample}_merged_combined.fastq 2> $2/${sample}_pychopper.log
    convertfasta2fastq $2/${sample}_merged_combined.fastq $2/${sample}_merged_combined.fasta
  
  fi
  
}

# 3) run_porechop <raw.fastq.gz> <output_dir>
# input: <raw>.fastq.gz
# output: <output_directory> with barcodes demultiplexed
run_porechop(){

    sample=$(basename $1)

    echo "Processing Sample $sample for Porechop"
    python ${PORECHOP} -i $1 -b $2 --format fastq --threads 16 \
      --check_reads 1000 \
      --discard_middle \
      --end_size 100 \
      --min_trim_size 15 \
      --extra_end_trim 1 \
      --end_threshold 75 \
      --verbosity 2
}


# 4) post_porechop_run_cutadapt <input_fastq> <output_dir>
post_porechop_run_cutadapt(){
  
  input_dir=$(dirname $1)
  name=$(basename $1 .fastq)
  
  if [ -f $2/${name}_combined.fastq ]; then
    echo $2/${name}_combined.fastq
    echo -e "Re-orientation and Cutadapt: \e[32mCompleted\e[0m"
  
  else

    # subset fastq file to polyA and polyT fasta (i.e. reads ending with PolyA and starting with polyT)
    # reads that end with AAAAAAAAAA = plus reads 
    # reads that start with TTTTTTTTTT = minus reads (need to be reverse complemented)
    echo "Subsetting fastq to polyA and polyT sequences"
    python ${SUBSETPOLYTAILS} --fa $1 --o_name ${name} --o_dir $2
    
    # working in output directory
    cd $2
    
    # reverse complement minus reads (reads ending with polyT)
    seqtk seq -r ${name}_PolyT.fastq > ${name}_PolyT_rev.fastq
    
    # use cutadapt package to trim polyA
    echo "Remove polyA sequences using cutadapt"
    cutadapt -a "A{60}" ${name}_PolyA.fastq -o ${name}_PolyA_cutadapted.fastq &> ${name}_polyA_cutadapt.log
    cutadapt -a "A{60}" ${name}_PolyT_rev.fastq -o ${name}_PolyT_rev_cuptadapted.fastq &> ${name}_polyT_cutadapt.log
    
    # concatenated reverse minus polyT and polyA reads
    cat ${name}_PolyA_cutadapted.fastq ${name}_PolyT_rev_cuptadapted.fastq > ${name}_combined.fastq
    
  fi
}


# 6) run_minimap2 <input_fasta> <output_dir>
# Aim: Align reads from trimming, filtering to genome of interest using Minimap2
# Input: <sample_name>_combined_reads.fasta
# Output: <sample_name>_combined_reads.sam, <sample_name>_Minimap2.log
run_minimap2(){

  name=$(basename $1 .fastq)

  if [ -f $2/${name}_sorted.sam ]; then
    echo -e "Minimap2: \e[32mCompleted\e[0m"
  
  else
  
    echo "Aligning ${name} using Minimap2"
    
    # remove secondary and supplementary alignments, but keep duplicates (required for phasing)
    minimap2 -t 46 -ax splice --secondary=no -R "@RG\tID:${name}\tSM:${name}\tLB:lib1\tPL:ONT" ${GENOME_FASTA} $1 > $2/${name}.sam 2> $2/${name}_minimap2.log
    samtools view -h -F 2308 $2/${name}.sam > $2/${name}_filtered.sam

    # sort sam file 
    samtools sort -O SAM $2/${name}_filtered.sam > $2/${name}_filtered_sorted.sam  

    # convert to bam file
    samtools view -S -b $2/${name}_filtered.sam | samtools sort -o $2/${name}_filtered_sorted.bam
    samtools index $2/${name}_filtered_sorted.bam
  
    htsbox samview -pS $2/${name}.sam > $2/${name}.paf
    awk -F'\t' '{if ($6!="*") {print $0}}' $2/${name}.paf > $2/${name}.filtered.paf
    awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $2/${name}.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $2/${name}"_mappedstats.txt"
  
  fi

}

run_minimap2stats(){

  name=$(basename $1 .fastq)

  htsbox samview -pS $2/${name}.sam > $2/${name}.paf
  awk -F'\t' '{if ($6!="*") {print $0}}' $2/${name}.paf > $2/${name}.filtered.paf
  awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $2/${name}.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $2/${name}"_mappedstats.txt"
  
}


# run_transcriptclean <input_sam> <output_dir>
run_transcriptclean(){
   
  name=$(basename $1 _merged_combined_filtered_sorted.sam)
  
  if [ -f $2/${name}/${name}_clean.TE.log ]; then
    echo -e "TranscriptClean: \e[32mCompleted\e[0m"
  
  else
  
    echo "TranscriptClean ${name}"  
    cd $2; mkdir -p ${name}
    cd $2/${name}
    python ${TCLEAN} --sam $1 --genome ${GENOME_FASTA} --outprefix $2/${name}/${name} --tmpDir $2/${name}/${name}_tmp --maxLenIndel=10
  
  fi
}


# 6) run_pbmm2 <input_fasta> <output_dir>
# Aim: re-align reads from transcript clean
# Input: <sample_name>_combined_reads.fasta
# Output: <sample_name>_combined_reads.sam, <sample_name>_Minimap2.log
run_pbmm2(){
  
  if [ -f $2/${name}_mapped.bam ]; then
    
    echo -e "Pbmm2: \e[32mCompleted\e[0m"
  
  else
  
    name=$(basename $1 _clean.fa)
    echo "TranscriptClean ${name}"
    echo "Aligning ${sample}: $1..."
    echo "Output: $2/${sample}_mapped.bam"
    
    cd $2
    pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} $1 ${name}_mapped.bam --log-level TRACE --log-file ${name}_mapped.log
    
  fi 
  
}


# filter_alignment <input_name> <input_mapped_dir>
filter_alignment(){

  if [ -f $2/$1.sorted.sam ]; then
  
    echo -e "Filtered: \e[32mCompleted\e[0m"
  
  else
    
    cd $2
    echo "Converting bam to sam and sort"
    samtools view -h $1.bam > $1.sam
    samtools bam2fq $1.bam| seqtk seq -A > $1.fa
    samtools sort -O SAM $1.sam > $1.sorted.sam
  
    # Alignment stats
    # Use the inforation in the paf file to create a new file where the columns correspond to the following: 
      #col1: name of the nanopore read 
      #col2: name of the sequence where nanopore read aligns (target sequence)
      #col3: start position of the alignment on the target sequence 
      #col4: length of the original nanopore read 
      #col5: length of the aligned part of the nanopore read  
      #col6: fraction of the aligned part of the nanopore read over the orginal length 
      #col7: fraction of the aligned part of the target sequence over the orginal length of the target sequence
      #col8: strand where the nanopore read aligns
      #col8: number of matched nucleotides of the nanopore read alignment on the target sequence
      #col9: identity (percentage of matched nucleotides over the aligned length of the nanopore read)
      #col10: number of mismatches of the nanopore read alignment on the target sequence
      #col11: number of insertions of the nanopore read alignment on the target sequence
      #col12: number of deletions of the nanopore read alignment on the target sequence
    
    echo "Dissecting alignment statistics"
    mkdir -p PAF; cd PAF
    htsbox samview -pS $2/$1.sorted.sam > $1.paf
    awk -F'\t' '{if ($6!="*") {print $0}}' $1.paf > $1.filtered.paf
    awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $1.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $1"_mappedstats.txt"
    ## filter based on alignable length (>0.85) and identity (>0.95)
    awk -F'\t' '{if ($6>=0.85 && $8>=0.95) {print $1}}' $1"_mappedstats.txt" > $1_filteredreads.txt
  
    picard FilterSamReads I=$2/$1.bam O=$2/$1.filtered.bam READ_LIST_FILE=$2/PAF/$1_filteredreads.txt FILTER=includeReadList &> $2/PAF/$1.picard.log
    samtools bam2fq $2/$1.filtered.bam| seqtk seq -A > $2/$1.filtered.fa
    samtools sort -O bam -o "$2/$1.filtered.sorted.bam" "$2/$1.filtered.bam"
    
    # https://bioinformatics.stackexchange.com/questions/3380/how-to-subset-a-bam-by-a-list-of-qnames
    #samtools view $2/$1.bam | grep -f $1_filteredreads.txt > $1.filtered.sam
    #samtools view -bS $1.filtered.sam > $1.filtered.bam
    #samtools bam2fq $2/$1.filtered.bam| seqtk seq -A > $2/$1.filtered.fa
  
  fi

}

# run_isoseq_collapse <input_aligned_bam> <output_name> <output_dir>
run_isoseq_collapse(){
      
  directory=$(dirname $1)
  
  if [ -f $directory/$2_collapsed.gff ]; then
    
    echo -e "Iso-Seq Collapse: \e[32mCompleted\e[0m"
  
  else
  
    echo "Collapsing..."
    echo "Output: $3/$2_collapsed.gff"
    
    cd ${directory}
    
    isoseq3 collapse $1 $2"_collapsed.gff" \
      --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
      --log-level TRACE --log-file $2"_collapsed.log"
  
  fi
}


# demuliplex_collapsed_isoforms <input_directory_fasta> <input_collapsed_directory> <output_name>
demuliplex_collapsed_isoforms(){
  
  if [ -f ${dir}/5_align/combined_fasta/$3"_sample_id.csv" ]; then
    
    echo -e "Extracted abundance: \e[32mCompleted\e[0m"
  
  else
  
    adapt_cupcake_to_ont.py $1 -o $3
  
    demux_cupcake_collapse.py \
      $2/$3"_collapsed.read_stat.txt" \
      ${dir}/5_align/combined_fasta/$3"_sample_id.csv"\
      --dataset=ont
  fi
  
}


# run_sqanti3 <gtf> <output_dir>
run_sqanti3(){
  
  if [ -f $2/${name}"_classification.txt" ]; then
  
    echo -e "SQANTI: \e[32mCompleted\e[0m"
  
  else
  
    name=$(basename $1 .gff)
  
    cd $2
   
    # sqanti qc
    echo "Processing Sample ${name} for SQANTI3 QC"
    python $SQANTI3_DIR/sqanti3_qc.py -v
    echo ${GENOME_GTF}
    echo ${GENOME_FASTA}
    
    python $SQANTI3_DIR/sqanti3_qc.py $1 ${GENOME_GTF} ${GENOME_FASTA} \
    --CAGE_peak ${CAGE_PEAK} \
    --polyA_motif_list ${POLYA} --skipORF \
    --genename --isoAnnotLite --report skip -t 30 &> ${name}.sqanti.qc.log
    
    echo "Processing Sample ${name} for SQANTI filter"
    python $SQANTI3_DIR/sqanti3_filter.py rules ${name}"_classification.txt" --gtf ${name}"_corrected.gtf" -j=${SQANTI_JSON} --skip_report &> ${name}.sqanti.filter.log
  
  fi
 
}

#### -------------------- post SQANTI -------------------

## ---------- run_cpat -----------------

# run_cpat <input_fasta> <output_name> <output_dir>
# Aim: 
  # call ORF from fasta file using CPAT (determine whether isoforms are protein-coding or non-protein-coding)
# Input:
  # input_fasta = input fasta for ORF to be called from
  # output_name = prefix output name
  # output_dir = path of output root directory to create CPAT folder directory
# Pre-requisite:
  # ${HEXAMER} = CPAT hexamer file (called from config file)
  # ${LOGITMODEL} = CPAT logit model (called from config file)
# Output
  # CPAT output files
  # CPAT log file

run_cpat(){
  
  mkdir -p $3/CPAT; cd $3/CPAT
  
  cpat.py --version
  cpat.py -x ${HEXAMER} -d ${LOGITMODEL} -g $1 --min-orf=50 --top-orf=50 -o $2 2> $2"_cpat.e"

  
}



## ---------- extract_best_orf -----------------

# extract_best_orf <output_name> <input/output_dir>
# Aim: 
  # extract the best ORF from CPAT for further analysis of ORF predictions for predicted NMD
# Input:
  # output_name = input and output name used from CPAT analysis
  # input/output_dir = directory of CPAT files
# Pre-requisite:
  # run_cpat to generate CPAT output files
# Output:
  # 

extract_best_orf(){
  
  cd $2/CPAT
  extract_fasta_bestorf.py --fa $1".ORF_seqs.fa" --orf $1".ORF_prob.best.tsv" --o_name $1"_bestORF" --o_dir $2 &> orfextract.log

}


## ---------- convert_gtf_bed12 -----------------

# convert_gtf_bed12 <input_gtf> 
# Aim:
  # convert gtf to bed12 file for downstream 
# Input: 
  # input_gtf = input gtf to be converted 
# Output:
  # path/to/original/directory/<sample>_sorted.bed12
convert_gtf_bed12(){
  
  # variables 
  output_dir="$(dirname $1)" 
  sample=${1%.gtf} # removes .gtf
  
  cd ${output_dir}
  
  gtfToGenePred $1 $sample.genePred
  genePredToBed $sample.genePred > $sample.bed12
  sort -k1,1 -k2,2n $sample.bed12 > $sample"_sorted.bed12"
  rm $sample.genePred $sample.bed12

}

## ---------- colour_by_abundance -----------------

# colour_by_abundance <cpat_name> <input_gtf> <input_abundance> <output_dir> <species=mouse/human>
# Aim:
  # generate multiple abundance file using input CPAT and expression using custom script 
  # custom script: ${LOGEN_ROOT}/merge_characterise_dataset/colour_transcripts_by_countandpotential.py
# Input:
  # cpat_name = CPAT input prefix names 
  # input_gtf = input gtf for conversion to bed12
  # input_abundance = path of input abundance file
  # output_dir = root of output directory for characterisation
  # species = mouse/human for determining CPAT threshold in script
# Output:
  # bed files
  
colour_by_abundance(){
  
  mkdir -p $4/bed12Files
  
  # convert gtf to bed12
  convert_gtf_bed12 $2
  
  # variables
  bed12=${2%.gtf}_sorted.bed12
  sample="$(basename $2)" 
  outputname=${sample%.gtf}
  echo $outputname
  
  colour_transcripts_by_countandpotential.py \
    --bed $bed12 \
    --cpat $4/CPAT/$1".ORF_prob.best.tsv" \
    --noORF $4/CPAT/$1".no_ORF.txt" \
    --a $3 \
    --o $outputname \
    --dir $4/bed12Files/ \
    --species $5
  
}


# subset_gene_reference <root_dir>
subset_gene_reference(){
  
  mkdir -p $1/TargetGenesRef
  
  source activate sqanti2_py3
  subset_reference_by_gene.py --r=${GENOME_GTF} --glist ${TGENES[@]} --o $1/TargetGenesRef
  
}


# run_transdecoder <name> <root_dir>
run_transdecoder(){
  source deactivate
  
  mkdir -p $2/Transdecoder; cd $2/Transdecoder
  
  TransDecoder.LongOrfs -t $2/CPAT/$1"_bestORF.fasta" &> transdecoder_longorf.log
  
  source activate sqanti2_py3
  hmmscan --cpu 8 --domtblout pfam.domtblout $PFAM_REF $1"_bestORF.fasta.transdecoder_dir"/longest_orfs.pep &> hmmscan.log
  
  source activate nanopore
  TransDecoder.Predict -t $2/CPAT/$1"_bestORF.fasta" --retain_pfam_hits pfam.domtblout --no_refine_starts &> transdecoder_predict.log
  sed '/^#/ d' < pfam.domtblout > pfam.domtblout.read
  
}