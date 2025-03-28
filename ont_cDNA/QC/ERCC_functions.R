
## ---------- packages -----------------

library("ggplot2")
library("dplyr")
library("stringr")
library("reshape2")
library("data.table")

#alignedstatsFiles <- "C:/Users/sl693/Dropbox/Scripts/LRPipeline/vignette/pilot/ERCC"
#ERCC_standard_dir <- "C:/Users/sl693/Dropbox/Scripts/LRPipeline/ont_cDNA/utils"
#SQANTI_dir <- "C:/Users/sl693/Dropbox/Scripts/LRPipeline/vignette/pilot/ERCC"


## ---------- functions -----------------

# read in mapped stats files from alignment to reference ERCCs
read_ERCC_alignedStats <- function(ERCC_alignedDir){
  colNames <- c(
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
  ERCC_alignedStats <- lapply(list.files(pattern = "mapped", path = ERCC_alignedDir, full.names = T), 
                              function(x) data.table::fread(x, col.names = colNames, data.table = F))
  names(ERCC_alignedStats) <- list.files(path = ERCC_alignedDir, pattern = "mapped")
  ERCC_alignedStats <- bind_rows(ERCC_alignedStats, .id = "file")
  ERCC_alignedStats$sample <- str_remove(ERCC_alignedStats$file, "_merged_combined_mappedstats.txt")
  
  return(ERCC_alignedStats)
}

# detect which ERCC is detected across dataset
detect_ERCC <- function(ERCC_alignedStats){
  # Create an empty list to store the results
  ercc_list <- lapply(unique(ERCC_alignedStats$sample), function(sample) {
    
    # Subset the `ercc_calculation` dataframe based on the current `sample`
    sample_ercc_calculation <- ercc_calculation
    
    # Create a new logical column indicating if the `ERCC_ID` is in the target sequences for this sample
    sample_ercc_calculation$barcode <- sample
    
    sample_ercc_calculation$Detected <- ifelse(sample_ercc_calculation$ERCC_ID %in% 
                                                 unique(ERCC_alignedStats[ERCC_alignedStats$sample == sample, "target_sequence"]), 
                                               TRUE, FALSE)
    
    # Return the modified dataframe for this sample
    return(sample_ercc_calculation)
  })
  
  ercc_combined <- bind_rows(ercc_list, .id = "sample") %>% 
    filter(!is.na(amount_of_ERCC)) %>%
    mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) 
  
  return(ercc_combined)

}

# plot the ERCC detected by amount and length
plot_amount_length <- function(ercc_combined) {
  
  create_plot <- function(y, y_label) {
    ggplot(ercc_combined, aes(x = Detected, y = !!sym(y), fill = Detected)) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.2)) +
      facet_grid(~sample) +
      theme_classic() +
      theme(
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none"
      ) +
      labs(x = NULL, y = y_label)
  }
  
  # Create the two plots
  p1 <- create_plot("log2_amount_of_ERCC", "Amount of ERCC (log2)")
  p2 <- create_plot("Length", "Length")
  
  return(list(p1, p2))
}

# plot the number of isoforms in classfile
plot_numIsoforms <- function(ERCC_classFile){
  p1 <- ERCC_classFile %>% group_by(chrom, structural_category) %>% tally() %>% 
    ggplot(., aes(x = chrom, y = n, fill = structural_category)) + 
    geom_bar(stat = "identity", position = position_dodge()) + 
    labs(x = "ERCC", y = "Number of unique isoforms") + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "bottom")
  
  p2 <- ERCC_classFile  %>% group_by(chrom) %>% tally() %>% 
    merge(., ercc_calculation, by.x = "chrom", by.y = "ERCC_ID", all = TRUE) %>% 
    mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) %>% 
    ggplot(., aes(x = n, y = log2_amount_of_ERCC)) + 
    geom_jitter(width = 0.2) + 
    labs(x = "Number of isoforms", y = "Amount of ERCC (Log2)") + 
    scale_y_continuous(trans='log10') +
    theme_classic()
  
  return(list(p1,p2))
}

plot_usage_persample <- function(ERCC_classFile, ERCC_demux) {
  
  # Merge data
  ERCC_classFileDemux <- merge(ERCC_classFile, ERCC_demux, by.x = "isoform", by.y = "id")
  
  plots <- list()
  percentages <- list()
  
  # Iterate over each sample column
  for (sample in colnames(ERCC_demux)[-1]) {
    dat <- ERCC_classFileDemux %>%
      select(all_of(c(sample, "isoform", "structural_category", "chrom"))) %>%
      filter(!!sym(sample) >= 1) %>%
      group_by(chrom) %>%
      mutate(
        perc = (!!sym(sample)) / sum(!!sym(sample)) * 100,
        major = ifelse(perc < 5, "minor", "major")
      ) %>%
      ungroup()
    
    # Separate major and minor isoforms
    major <- dat %>% filter(major != "minor") %>% select(chrom, perc, isoform, structural_category)
    minor <- dat %>% filter(major == "minor")
    
    # Group minor isoforms and sum percentages
    minorgrouped <- aggregate(minor$perc, by = list(chrom = minor$chrom), FUN = sum) %>%
      mutate(isoform = "minor", structural_category = "minor") %>%
      dplyr::rename(perc = x)
    
    # Tally minor isoforms
    minortally <- minor %>% group_by(chrom) %>% tally() %>% mutate(isoform = "minor")
    
    # Combine data for table output
    percentage_table <- rbind(major, minorgrouped) %>%
      full_join(minortally, by = c("isoform", "chrom")) %>%
      arrange(chrom, desc(perc))
    
    # Store the table
    percentages[[sample]] <- percentage_table
    
    # Create the plot
    percentage_table$structural_category <- factor(percentage_table$structural_category, 
                                               levels = rev(unique(percentage_table$structural_category)))

    plots[[sample]] <- percentage_table %>%
      ggplot(aes(x = chrom, y = as.numeric(perc), fill = structural_category)) +
      geom_bar(stat = "identity", color = "black", size = 0.2) +
      theme_classic() +
      labs(x = "Gene", y = "Isoform fraction (%)", title = sample) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(
        name = "Isoform Classification",
        values = rev(c(alpha("#00BFC4", 0.8), alpha("#00BFC4", 0.3),
                       alpha("#F8766D", 0.8), alpha("#F8766D", 0.3),
                       alpha("#808080", 0.3)))
      ) +
      geom_text(aes(label = n), color = "black", size = 4, position = position_stack(vjust = 0.5)) +
      theme(legend.position = "bottom")
  }
  
  # Return both plots and tables
  return(list(plots = plots, tables = percentages))
}

tabulate_isoformFraction <- function(pUsageOutput){
  
  # list the output from isoform fraction function and remove minor isoforms
  dat <- rbindlist(pUsageOutput$tables, idcol = "sample") %>% 
    select(sample, chrom, perc, isoform) %>% 
    filter(isoform != "minor")
  
  # tally the number of major isofoms per sample per ERCC
  dat2 <- dat %>% group_by(sample, chrom) %>% tally()
  
  # table of ERCC by sample and fraction 
  out <- merge(dat,dat2, by = c("sample","chrom"))
  colnames(out) <- c("sample","ERCC","Isoform fraction (%)", "Isoform", "Number of isoforms")
  
  return(out)
}

