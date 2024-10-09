
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

  source activate nanopore
  seqtk seq -A $1 > $2
}

merge_fastq_across_samples(){
  # variables 
  gval=$1
  input_dir=$2
  output_dir=$3
  
  echo "Merging ${gval}"
  
  fastq=$(ls ${input_dir}/${gval}/*fa* 2>/dev/null)
  num_files=$(echo "$fastq" | wc -w)
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

  source activate nanopore
  seqkit stats ${output_dir}/${gval}_merged.fastq > ${output_dir}/${gval}_readstats.txt
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
  
  source activate nanopore 
  
  # requires fasta files for downstream
  echo "Converting $1 to fasta"
  seqtk seq -a $1 > ${input_dir}/${name}.fasta
  
  # subset fasta file to polyA and polyT fasta (i.e. reads ending with PolyA and starting with polyT)
  # reads that end with AAAAAAAAAA = plus reads 
  # reads that start with TTTTTTTTTT = minus reads (need to be reverse complemented)
  echo "Subsetting fasta to polyA and polyT sequences"
  python ${SUBSETPOLYTAILS} --fa ${input_dir}/${name}.fasta --o_name ${name} --o_dir $2
  
  # working in output directory
  cd $2
  
  # reverse complement minus reads (reads ending with polyT)
  seqtk seq -r ${name}_PolyT.fasta > ${name}_PolyT_rev.fasta
  
  # use cutadapt package to trim polyA
  echo "Remove polyA sequences using cutadapt"
  cutadapt -a "A{60}" ${name}_PolyA.fasta -o ${name}_PolyA_cutadapted.fasta &> ${name}_polyA_cutadapt.log
  cutadapt -a "A{60}" ${name}_PolyT_rev.fasta -o ${name}_PolyT_rev_cuptadapted.fasta &> ${name}_polyT_cutadapt.log
  
  # concatenated reverse minus polyT and polyA reads
  cat ${name}_PolyA_cutadapted.fasta ${name}_PolyT_rev_cuptadapted.fasta > ${name}_combined.fasta
  
  source deactivate
}


# 6) run_minimap2 <input_fasta> <output_dir>
# Aim: Align reads from trimming, filtering to genome of interest using Minimap2
# Input: <sample_name>_combined_reads.fasta
# Output: <sample_name>_combined_reads.sam, <sample_name>_Minimap2.log
run_minimap2(){

  source activate nanopore
  name=$(basename $1 .fasta)
  echo "Aligning ${name} using Minimap2"
  
  minimap2 -t 46 -ax splice ${GENOME_FASTA} $1 > $2/${name}.sam 2> $2/${name}_minimap2.log
  samtools sort -O SAM $2/${name}.sam > $2/${name}_sorted.sam

  htsbox samview -pS $2/${name}_sorted.sam > $2/${name}.paf
  awk -F'\t' '{if ($6!="*") {print $0}}' $2/${name}.paf > $2/${name}.filtered.paf
  awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $2/${name}.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $2/${name}"_mappedstats.txt"

}

# run_transcriptclean <input_sam> <output_dir>
run_transcriptclean(){
  
  source activate sqanti2_py3
   
  name=$(basename $1 _merged_combined_sorted.sam)
  echo "TranscriptClean ${name}"
  
  cd $2; mkdir -p ${name}
  cd $2/${name}
  python ${TCLEAN} --sam $1 --genome ${GENOME_FASTA} --outprefix $2/${name}/${name} --tmpDir $2/${name}/${name}_tmp
}


# 6) run_pbmm2 <input_fasta> <output_dir>
# Aim: re-align reads from transcript clean
# Input: <sample_name>_combined_reads.fasta
# Output: <sample_name>_combined_reads.sam, <sample_name>_Minimap2.log
run_pbmm2(){
  
  source activate isoseq3
  name=$(basename $1 _clean.fa)
  echo "TranscriptClean ${name}"
  echo "Aligning ${sample}: $1..."
  echo "Output: $2/${sample}_mapped.bam"
  
  cd $2
  pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} $1 ${name}_mapped.bam --log-level TRACE --log-file ${name}_mapped.log
  
}


# filter_alignment <input_name> <input_mapped_dir>
filter_alignment(){
  
  source activate nanopore
  
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

  source activate sqanti2
  picard FilterSamReads I=$2/$1.bam O=$2/$1.filtered.bam READ_LIST_FILE=$2/PAF/$1_filteredreads.txt FILTER=includeReadList &> $2/PAF/$1.picard.log
  
  source activate nanopore
  samtools bam2fq $2/$1.filtered.bam| seqtk seq -A > $2/$1.filtered.fa
  samtools sort -O bam -o "$2/$1.filtered.sorted.bam" "$2/$1.filtered.bam"
  
  # https://bioinformatics.stackexchange.com/questions/3380/how-to-subset-a-bam-by-a-list-of-qnames
  #source activate nanopore
  #samtools view $2/$1.bam | grep -f $1_filteredreads.txt > $1.filtered.sam
  #samtools view -bS $1.filtered.sam > $1.filtered.bam
  #samtools bam2fq $2/$1.filtered.bam| seqtk seq -A > $2/$1.filtered.fa

}

# run_isoseq_collapse <input_aligned_bam> <output_name> <output_dir>
run_isoseq_collapse(){
  echo "Collapsing..."
  echo "Output: $3/$2_collapsed.gff"
  
  directory=$(dirname $1)
  cd ${directory}

  source activate isoseq3
  
  isoseq3 collapse $1 $2"_collapsed.gff" \
    --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
    --log-level TRACE --log-file $2"_collapsed.log"
}


# demuliplex_collapsed_isoforms <input_directory_fasta> <input_collapsed_directory> <output_name>
demuliplex_collapsed_isoforms(){
  adapt_cupcake_to_ont.py $1 -o $3

  demux_cupcake_collapse.py \
    $2/$3"_collapsed.read_stat.txt" \
    ${dir}/5_align/combined_fasta/$3"_sample_id.csv"\
    --dataset=ont
  
}


# run_sqanti3 <gtf> <output_dir>
run_sqanti3(){
  
  name=$(basename $1 .gff)

  cd $2
  source activate sqanti2_py3
  
  # sqanti qc
  echo "Processing Sample ${name} for SQANTI3 QC"
  python $SQANTI3_DIR/sqanti3_qc.py -v
  echo ${GENOME_GTF}
  echo ${GENOME_FASTA}
  
  python $SQANTI3_DIR/sqanti3_qc.py $1 ${GENOME_GTF} ${GENOME_FASTA} \
  --CAGE_peak ${CAGE_PEAK} \
  --polyA_motif_list ${POLYA} \
  --genename --isoAnnotLite --report skip -t 30 &> ${name}.sqanti.qc.log
  
  echo "Processing Sample ${name} for SQANTI filter"
  python $SQANTI3_DIR/sqanti3_filter.py rules ${name}"_classification.txt" ${name}"_corrected.fasta" ${name}"_corrected.gtf" -j=${filteringJson} --skip_report &> ${name}.sqanti.filter.log

  
}
