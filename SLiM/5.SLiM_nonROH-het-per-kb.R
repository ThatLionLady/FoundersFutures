library(data.table)
library(dplyr)
library(ggplot2)

# Uses: per cycle roh file ex. ROH/ROH_A1_N0_10.hom
#       per sample het output ex. VCFtools/noROH_indiv_A1_N0_10_i0.het
#       per sample bed file ex. BED/roh_A1_N0_10_i0.bed

setwd("ZDC/SLiM")

# Define chromosomes and simulation directories
chromosomes <- c("A1")
sim_dirs <- c("N", "W", "R", "D")

# Load chromosome lengths from reference file
fai <- fread('References/LionRefs/LionRef.fasta.fai', header=FALSE)
colnames(fai) <- c('Chr', 'Chrlgth', 'Offset', 'Base', 'Width')

# Initialize an empty list to store all results
het_results_list <- list()

# Loop through chromosomes, simulations, and cycles
for (chr in chromosomes) {
  cat(sprintf("Processing Chromosome: %s at %s\n", chr, Sys.time()))
  chr_length <- fai$Chrlgth[fai$Chr == chr]
  genome_length_kb <- chr_length / 1000  # Convert chromosome length to kilobases
  
  for (sim_dir in sim_dirs) {
    sim_name <- toupper(sim_dir)
    cat(sprintf("  Simulation: %s\n", sim_name))
    
    for (run in 0:9) { # Assuming 10 runs
      cat(sprintf("    Processing Run %d for Chromosome %s, Simulation %s\n", run, chr, sim_name))

      for (cycle in seq(5, 20, 5)) { # feature every 5 cycles for 20 cycles
        # Define the path to the ROH file
        sample_file <- paste0("Lists/samples_", chr, "_", sim_name, run, "_", cycle, ".txt")
        
        # Check if the file exists
        if (file.exists(sample_file)) {
          cat(sprintf("      Processing Run %d for Chromosome %s, Simulation %s Cycle %d\n", run, chr, sim_name, cycle))
          
          # Loop through samples and compute heterozygosity per kb
          samples <- readLines(sample_file)
          total_samples <- length(samples)
          sample_counter <- 1

          for (SAMPLE in samples) {
            cat(sprintf("        Processing Sample %d of %d: %s\n", sample_counter, total_samples, SAMPLE))
            sample_counter <- sample_counter + 1
            
            het_file <- paste0("HET/noROH_indiv_", chr, "_", sim_name, run, "_", cycle, "_", SAMPLE, ".het")
            roh_bed_file <- paste0("BED/roh_", chr, "_", sim_name, run, "_", cycle,  "_", SAMPLE, ".bed")
            
            # Check if het file exists and is non-empty
            if (!file.exists(het_file) || file.info(het_file)$size == 0) {
              cat(sprintf("          Skipping Sample %s: het file missing or empty\n", SAMPLE))
              next
            }
            
            het_df <- read.table(het_file, header = TRUE)
            het_df <- het_df[het_df$INDV == SAMPLE, ]
            if (nrow(het_df) == 0) {
              cat(sprintf("          Skipping Sample %s: no matching INDV\n", SAMPLE))
              next
            }
            
            # Default to full genome length
            non_roh_length_kb <- genome_length_kb
            
            # Try to read the BED file
            if (file.exists(roh_bed_file)) {
              roh_bed_df <- tryCatch(
                {
                  read.table(roh_bed_file, header = FALSE, sep = "\t")
                },
                error = function(e) {
                  cat(sprintf("          Error reading BED file for sample %s: %s\n", SAMPLE, e$message))
                  return(NULL)
                }
              )
              
              if (!is.null(roh_bed_df) && nrow(roh_bed_df) > 0) {
                non_roh_length_kb <- genome_length_kb - sum(roh_bed_df$V3 - roh_bed_df$V2) / 1000
              } else {
                cat(sprintf("          No ROH entries found for sample %s — using full genome length\n", SAMPLE))
              }
            } else {
              cat(sprintf("          BED file missing for sample %s — using full genome length\n", SAMPLE))
            }
            
            # Calculate heterozygosity per kb
            het_df$HET_PER_KB <- (het_df$N_SITES - het_df$O.HOM.) / non_roh_length_kb
            het_df$non_roh_length_kb <- non_roh_length_kb
            
            # Add metadata
            het_df$Sample <- SAMPLE
            het_df$Simulation <- sim_name
            het_df$Run <- run
            het_df$Cycle <- cycle
            het_df$Chromosome <- chr
            
            # Append to results
            het_results_list <- append(het_results_list, list(het_df))
          } # end sample loop
        } else {
          cat(sprintf("      Missing files for Cycle %d, Run %d: %s or %s\n", cycle, run, het_file, roh_bed_file))
       } # end error from file check
      } # end cycle loop
    } # end run lopp
  } # end simulation loop
} # end chromosome loop

# Combine all results into a single data.table
het_results <- rbindlist(het_results_list)

# Save all results to a single CSV file
output_file <- "heterozygosity_per_kb_all_samples.csv"
write.table(het_results, file=output_file, sep=",", quote=FALSE, row.names=FALSE)

het_results <- het_results %>%
  mutate(Simulation = factor(Simulation, levels = c("R", "N", "W", "D")))

# Check for Variation
het_results %>%
  group_by(Simulation, Run, Cycle) %>%
  summarise(sd_het = sd(HET_PER_KB, na.rm = TRUE), .groups = "drop") %>%
  summarise(min_sd = min(sd_het), max_sd = max(sd_het))

# Compute per-run means
het_run_means <- het_results %>%
  group_by(Simulation, Cycle, Run) %>%
  summarise(
    Run_Mean_HET = mean(HET_PER_KB, na.rm = TRUE),
    Run_SD_HET = sd(HET_PER_KB, na.rm = TRUE),
    .groups = "drop"
  )

# Between-run variation
het_between <- het_run_means %>%
  group_by(Simulation, Cycle) %>%
  summarise(
    Mean_HET = mean(Run_Mean_HET, na.rm = TRUE),
    SD_between = sd(Run_Mean_HET, na.rm = TRUE),
    n_runs = n(),
    SE = SD_between / sqrt(n_runs),
    CI_Lower = Mean_HET - qt(0.975, df = n_runs - 1) * SE,
    CI_Upper = Mean_HET + qt(0.975, df = n_runs - 1) * SE,
    .groups = "drop"
  )

# Plot
ggplot() +
  
  # each run trajectory
  geom_line(
    data = het_run_means,
    aes(x = Cycle, y = Run_Mean_HET, group = Run, color = Simulation),
    alpha = 0.25
  ) +
  
  # mean across runs
  geom_line(
    data = het_between,
    aes(x = Cycle, y = Mean_HET, color = Simulation),
    linewidth = 1.2
  ) +
  
  # between-run CI
  geom_ribbon(
    data = het_between,
    aes(x = Cycle, ymin = CI_Lower, ymax = CI_Upper, fill = Simulation),
    alpha = 0.2,
    color = NA
  ) +
  
  facet_wrap(~ Simulation) +
  theme_minimal() +
  labs(
    title = "Between-Run Variation in Heterozygosity (HET_PER_KB)",
    x = "Cycle",
    y = "Heterozygosity per kb"
  )

ggsave("Plots/HET_PER_KB_scale_Run_variation.pdf", width = 8, height = 6)
ggsave("Plots/HET_PER_KB_scale_Run_variation.png", width = 8, height = 6)
