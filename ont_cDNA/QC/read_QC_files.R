#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Purpose: read files for QC report 
## Output: QC_Input.RData in --rootdir
## Rscript read_QC_files.R --m manifest.csv --r root_dir
## --------------------------------


## ---------- packages -----------------

# List of required packages
required_packages <- c("data.table", "dplyr", "stringr", "reshape2", "optparse")

# Check if packages are installed, and install missing ones
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    suppressMessages(install.packages(pkg, repos = 'https://cran.uk.r-project.org'))
    suppressWarnings(install.packages(pkg, repos = 'https://cran.uk.r-project.org'))
  }
}

# load packages
suppressMessages({
  lapply(required_packages, library, character.only = TRUE)
})


## ---------- arguments -----------------

# Define the option list
option_list <- list( 
  make_option(c("-m", "--manifest"), type="character", help="manifest file", metavar="character"),
  make_option(c("-r", "--rootdir"), type="character", help="root working directory containing files", metavar="character")
)

# Parse the command-line options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate that required arguments are provided
if (is.null(opt$manifest)) {
  stop("Error: '--manifest' argument is required. Use -m or --manifest to specify the manifest file.")
}

if (is.null(opt$rootdir)) {
  stop("Error: '--rootdir' argument is required. Use -r or --rootdir to specify the root directory.")
}

# Print the parsed arguments 
message("Manifest file: ", opt$manifest)
manifest <- read.csv(opt$manifest, stringsAsFactors = FALSE)
message("Root directory: ", opt$rootdir)
rootDir <- opt$rootdir

## ---------- functions -----------------

# check if correctly labelled according to manifest 
check_sample <- function(Files){
  ID <- word(basename(Files),c(1), sep = fixed("_"))
  if(length(intersect(manifest$sample, ID)) == 0){
    message("Error: no matching samples between manifest and files")
    message("manifest sample column:", paste0(manifest$sample, sep = ","))
    message("file names:", paste0(ID , sep = ","))
    stop("Exiting script due to mismatched samples between alignment files and manifest.")
  }
}

check_id <- function(demux){
  if(length(intersect(manifest$ID, colnames(demux))) == 0){
    message("Error: no matching IDs between manifest and demux file column")
    message("manifest id column:", paste0(manifest$ID, sep = ","))
    message("column names of demux count file:", paste0(colnames(demux), sep = ","))
    stop("Exiting script due to mismatched samples between demux expression file and manifest.")
  }
}

## ---------- manifest file -----------------
# manifest individual as character
manifest$individual <- as.character(manifest$individual)

## ---------- demultiplex stats -----------------

message("Reading in demultipled read stats")
demuxFiles <- list.files(path = paste0(rootDir, "/1b_demultiplex_merged/"), pattern = "readstats", full = T, recursive = F)
print(demuxFiles)
check_sample(demuxFiles)

# read in files
demuxStats <- lapply(demuxFiles, function(x) fread(x, data.table = F))
names(demuxStats) <- list.files(path = paste0(rootDir, "/1b_demultiplex_merged/"), pattern = "readstats", full = F)
demuxStats <- rbindlist(demuxStats, idcol = "sampleNum") %>% mutate(sample = word(sampleNum,c(1),sep=fixed("_")))
demuxStats$num_seqs <- as.numeric(gsub(",", "", demuxStats$num_seqs))
demuxStats <- merge(demuxStats, manifest, by = "sample") 

## ---------- alignment stats -----------------

# Column names for alignment stats
nameAlignedStats <- c(
  "read_name",                # Name of the nanopore read
  "target_sequence",          # Name of the target sequence
  "start_position",           # Start position of the alignment on the target sequence
  "read_length",              # Length of the original nanopore read
  "aligned_length",           # Length of the aligned part of the nanopore read
  "aligned_fraction",         # Fraction of the aligned part of the nanopore read over the original length
  "matches",                  # Number of matched nucleotides in the alignment
  "identity_percentage",      # Fraction of matched nucleotides over the aligned length
  "strand",                   # Strand where the nanopore read aligns
  "mismatches",               # Number of mismatches in the alignment
  "insertions",               # Number of insertions in the alignment
  "deletions"                 # Number of deletions in the alignment
)


read_alignment_stats <- function(path, pattern, col_names, id_column_name = "file") {
  message("Reading alignment stats from: ", path)
  
  # List files matching the pattern
  file_list <- list.files(path = path, pattern = pattern, full.names = TRUE, recursive = FALSE)
  
  # print the files read
  print(file_list)
  
  # Check file IDs
  check_sample(file_list)
  
  # Read and process alignment stats
  alignment_data <- lapply(file_list, function(x) fread(x, col.names = col_names, data.table = FALSE))
  
  # Assign file names as IDs
  names(alignment_data) <- list.files(path = path, pattern = pattern, full.names = FALSE)
  
  # Combine into a single data table
  alignment_data <- rbindlist(alignment_data, idcol = id_column_name)
  
  # Extract sample name
  alignment_data$sample <- sub("_.*", "", alignment_data[[id_column_name]])
  
  return(alignment_data)
}

# Alignment stats from Minimap
mappedStats <- read_alignment_stats(
  path = paste0(rootDir, "/3_minimap/"), 
  pattern = "mapped", 
  col_names = nameAlignedStats
)

# Alignment stats from pbmm2
filteredAlignment <- read_alignment_stats(
  path = paste0(rootDir, "/5_cupcake/5_align/PAF/"), 
  pattern = "mappedstats", 
  col_names = nameAlignedStats
)

## ---------- SQANTI classification file -----------------

# classifiction file with all isoforms and artifacts
preFilteredClassFileName <- list.files(path = paste0(rootDir,"/5_cupcake/7_sqanti3"), pattern = "RulesFilter_result_classification.txt$", full = T)
message("Reading SQANTI classification file: ", preFilteredClassFileName)
preFilteredClassFile <- fread(preFilteredClassFileName, data.table = F)

# datawrangle columns
xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog","genic","antisense","fusion","intergenic","genic_intron")
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron")
preFilteredClassFile$structural_category = factor(preFilteredClassFile$structural_category,
                                                  labels = xaxislabelsF1, 
                                                  levels = xaxislevelsF1,
                                                  ordered=TRUE)

# retain only isoforms kept from SQANTI filtering
class_file <- preFilteredClassFile %>% filter(filter_result == "Isoform")

# filter mono-exonic isoforms
monoexonic <- class_file[class_file$subcategory == "mono-exon" & class_file$structural_category != "FSM","isoform"]
class_file <- class_file %>% filter(!isoform %in% monoexonic)


## ---------- demux expression file -----------------

demux <- list.files(path = paste0(rootDir,"/5_cupcake/6_collapse"), pattern = "count.csv", full = T)
message("Reading demux expression file: ", demux)
demux <- fread(demux, data.table = F)
check_id(demux)
demuxFiltered  <- demux %>% filter(id %in% class_file$isoform) %>% tibble::column_to_rownames(., var = "id")

## ---------- SQANTI filtering reasons file -----------------

ReasonFiltered <- list.files(path = paste0(rootDir,"/5_cupcake/7_sqanti3"), pattern = "filtering_reasons", full = T)
message("Reading SQANTI Reasons filtered file: ", ReasonFiltered)
ReasonFiltered  <- fread(ReasonFiltered, data.table = F)

## ---------- Output to Rdata -----------------

input <- list(
  manifest = manifest, 
  demuxStats = demuxStats,
  mappedStats = mappedStats, 
  filteredAlignment = filteredAlignment,
  preFilteredClassFile = preFilteredClassFile,
  class_file = class_file,
  demux = demux, 
  demuxFiltered = demuxFiltered,
  ReasonFiltered = ReasonFiltered
)

message("Write output to: ", paste0(rootDir, "/QC_input.RData"))
save(input, file = paste0(rootDir, "/QC_input.RData"))
