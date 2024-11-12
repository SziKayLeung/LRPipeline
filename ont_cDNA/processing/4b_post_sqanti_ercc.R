library("ggplot2")
library("dplyr")
library("stringr")
library("reshape2")

demuxStats <- lapply(list.files(pattern = "demux", full.names = T), function(x) read.table(x, header = T))
demuxStats <- do.call("rbind", demuxStats) %>% dplyr::mutate(barcode = stringr::word(file, c(1), sep = fixed("_")),
                                                             stage = "demuxed")

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
alignedStats <- lapply(list.files(pattern = "mapped"), function(x) read.table(x, col.names = nameAlignedStats))
alignedStats <- bind_rows(alignedStats, .id = "barcode")

demuxStats[,c("barcode","stage","num_seqs")]

detected <- list(
  barcode1 = unique(alignedStats[alignedStats$barcode == 1,"target_sequence"]),
  barcode2 = unique(alignedStats[alignedStats$barcode == 2,"target_sequence"]),
  barcode3 = unique(alignedStats[alignedStats$barcode == 3,"target_sequence"])
)

ercc_calculation <- read.csv("ERCC_calculations.csv")
ercc_stats <- read.csv("ERCC_Stats.csv")

ercc_calculation_df <- ercc_calculation %>% mutate(barcode1 = ifelse(ERCC_ID %in% detected$barcode1, TRUE, FALSE),
                            barcode2 = ifelse(ERCC_ID %in% detected$barcode2, TRUE, FALSE),
                            barcode3 = ifelse(ERCC_ID %in% detected$barcode3, TRUE, FALSE)) %>% 
  select(ERCC_ID, amount_of_ERCC, barcode1, barcode2, barcode3) %>% 
  reshape2::melt(id = c("ERCC_ID", "amount_of_ERCC"), value.name = "Detected", variable.name = "barcode") %>% 
  filter(!is.na(amount_of_ERCC)) %>%
  mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) 
ggplot(., aes(x = Detected, y = log2_amount_of_ERCC)) + geom_boxplot() + 
  geom_point(aes(colour = ERCC_ID), position = position_jitter(width = 0.2), alpha = 0.5) +
  facet_grid(~barcode) + theme_classic()

ercc_stats_df <- ercc_stats  %>% mutate(barcode1 = ifelse(ERCC.ID %in% detected$barcode1, TRUE, FALSE),
                            barcode2 = ifelse(ERCC.ID %in% detected$barcode2, TRUE, FALSE),
                            barcode3 = ifelse(ERCC.ID %in% detected$barcode3, TRUE, FALSE)) %>% 
  select(ERCC.ID, Length, barcode1, barcode2, barcode3) %>% 
  reshape2::melt(id = c("ERCC.ID", "Length"), value.name = "Detected", variable.name = "barcode") %>%
  mutate(ERCC_ID = ERCC.ID)

  ggplot(., aes(x = Detected, y = Length)) + geom_boxplot() + 
  geom_point(aes(colour = ERCC.ID), position = position_jitter(width = 0.2), alpha = 0.5) +
  facet_grid(~barcode) + theme_classic()
  
merge(ercc_stats_df, ercc_calculation_df, by = c("ERCC_ID", "barcode")) %>% 
  mutate(Detected.x = factor(Detected.x, levels = c("TRUE","FALSE"))) %>%
  ggplot(aes(y = log2_amount_of_ERCC, x = Length, colour = Detected.x)) + geom_point() + facet_grid(~barcode) +
  labs(x = "known length (bp)", y = "Amount of ERCC (Log2)", colour = "Detected") +
  theme_classic()



alignedStats %>% group_by(barcode) %>% tally()
ggplot(alignedStats, aes(x = aligned_fraction, y = identity_percentage, colour = target_sequence)) + 
  geom_point() +
  xlim(0,1) + 
  facet_grid(~barcode) +
  geom_vline(xintercept=0.85, linetype = "dotted") +
  geom_hline(yintercept=0.95, linetype = "dotted") +
  theme_classic() +
  labs(x = "Alignment Coverage", y = "Alignment Identity")

# ERCC
cat("Total unique ERCCs:", length(unique(alignedStats$target_sequence)), paste0("(", round(length(unique(alignedStats$target_sequence))/92 * 100,2), "%)"))

# redundant 
redundant <- alignedStats %>% group_by(barcode, target_sequence) %>% tally()
colnames(redundant) <- c("barcode","ERCC","num_isoforms")

# isoform vs concentration
isoform_conc <- merge(redundant, ercc_calculation, by.x = "ERCC", by.y = "ERCC_ID", all = TRUE) %>%
  mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) %>%
  tidyr::replace_na(list(num_isoforms = 0)) %>% 
  filter(!is.na(amount_of_ERCC), !is.na(barcode))

ggplot(isoform_conc, aes(x = num_isoforms, y = log2_amount_of_ERCC, colour = barcode)) + 
  geom_jitter(width = 0.2) + 
  labs(x = "Number of reads", y = "Amount of ERCC (Log2)") + 
  scale_y_continuous(trans='log10') +
  theme_classic()

collapsed <- read.table("sfariReview_collapsed.read_stat.txt", header = T)

merge(collapsed, alignedStats, by.x = "id", by.y = "read_name", all.x = T) 

alignedStats[alignedStats$aligned_fraction  >= 0.85,]
nrow(alignedStats)
