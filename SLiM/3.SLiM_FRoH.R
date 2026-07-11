library(data.table)
library(dplyr)
library(ggplot2)

setwd("ZDC/SLiM/")

# Define chromosomes and simulation directories
chromosomes <- c("A1")
sims <- c("N", "W", "R", "D")

# Load chromosome lengths from reference file
fai <- fread('References/LionRefs/LionRef.fasta.fai', header=FALSE)
colnames(fai) <- c('Chr', 'Chrlgth', 'Offset', 'Base', 'Width')

# Initialize final results data frame
#final_results <- data.frame(Simulation=character(), Chromosome=character(), Cycle=integer(), Mean_FROH=numeric(), CI_Lower=numeric(), CI_Upper=numeric(), stringsAsFactors=FALSE)
final_results <- data.frame(Sample=character(),
                            Chromosome=character(),
                            Simulation=character(),
                            Cycle=integer(),
                            Run=integer(),
                            FROH_2Mb=numeric(),
                            FROH_1Mb=numeric(),
                            FROH_50Kb=numeric(),
                            FROH_all=numeric(),
                            stringsAsFactors=FALSE)

# Create an empty combined dataframe to accumulate ROH data
roh_combined <- data.frame(FID=character(), ID=character(), PHE=numeric(), Chr=character(), SNP1=numeric(), SNP2=numeric(), Start=numeric(), End=numeric(), KB=numeric(), NSNP=numeric(), DENSITY=numeric(), PHOM=numeric(), PHET=numeric(), Length=numeric(), Length_Mb=numeric(), Category=character(), stringsAsFactors=FALSE)

# Loop through chromosomes, simulations, and cycles
for (chr in chromosomes) {
  cat(sprintf("Processing Chromosome: %s at %s\n", chr, Sys.time()))
  chr_length <- fai$Chrlgth[fai$Chr == chr]
  
  for (sim_name in sims) {
    cat(sprintf("  Simulation: %s\n", sim_name))
    
    froh_data <- data.frame()
    
    for (run in 0:9) { # 10 runs labeled 0-9
      cat(sprintf("    Processing Run %d for Chromosome %s, Simulation %s\n", run, chr, sim_name))
      
      for (cycle in seq(5, 20, 5)) { # feature every 5 cycles for 20 cycles
        # Define the path to the ROH file
        roh_file <- paste0("ROH/ROH_", chr, "_", sim_name, run, "_", cycle, ".hom")
        
        # Check if the file exists
        if (file.exists(roh_file)) {
          cat(sprintf("      Processing Cycle %d for Run %d\n", cycle, run))
          
          # Read the ROH data
          roh <- fread(roh_file, header=TRUE)
          colnames(roh) <- c('FID', 'ID', 'PHE', 'Chr', 'SNP1', 'SNP2', 'Start', 'End', 'KB', "NSNP", "DENSITY", "PHOM", "PHET")
          roh$Length <- roh$KB * 1000 # Convert Kb to bases
          
          # Calculate FROH for each sample
          sample_list <- unique(roh$ID)
          total_samples <- length(sample_list)
          sample_counter <- 1

          for (SAMPLE in sample_list) {
            cat(sprintf("        Processing Sample %d of %d: %s\n", sample_counter, total_samples, SAMPLE))
            sample_counter <- sample_counter + 1
            
            roh_sample <- roh[roh$ID == SAMPLE, ]
            
            # Calculate FROH > 2Mb
            froh_2Mb <- sum(roh_sample$Length[roh_sample$Length > 2000000], na.rm=TRUE) / chr_length

            # Calculate FROH > 1Mb
            froh_1Mb <- sum(roh_sample$Length[roh_sample$Length > 1000000], na.rm=TRUE) / chr_length
            
            # Calculate FROH > 50Kb
            froh_50Kb <- sum(roh_sample$Length[roh_sample$Length > 50000], na.rm=TRUE) / chr_length
            
            # Calculate FROH for all ROHs
            froh_all <- sum(roh_sample$Length, na.rm=TRUE) / chr_length
            
            cat(sprintf("Sample: %s, Chromosome: %s, Simulation: %s, Cycle: %d, Run: %d\n", 
                        SAMPLE, chr, sim_name, cycle, run))
            cat(sprintf("FROH_2Mb: %.7f, FROH_50Kb: %.7f, FROH_all: %.7f\n", 
                        froh_2Mb, froh_50Kb, froh_all))
            
            # Append to froh_data
            froh_data <- rbind(froh_data, data.frame(
              Sample = SAMPLE,
              Chromosome = chr,
              Simulation = sim_name,
              Cycle = cycle,
              Run = run,
              FROH_2Mb = as.numeric(froh_2Mb),
              FROH_1Mb = as.numeric(froh_1Mb),
              FROH_50Kb = as.numeric(froh_50Kb),
              FROH_all = as.numeric(froh_all)))

            ################################################ Combine Data
      
            # Combine the current sample's data with the combined dataframe
            roh_combined <- rbind(roh_combined, roh_sample)
            
            ################################################ Create BED File

            roh_sample <- roh_sample[Length > 2000000]
            
            # Convert to BED format
            setDT(roh_sample)
            bed_df <- roh_sample[, .(Chr, Start, End)]
            setnames(bed_df, c("Chr", "Start", "End"), c("#CHR", "START", "END"))
            
            # Save BED file
            roh_bed_file <- paste0("BED/roh_", chr, "_", sim_name, run, "_", cycle, "_", SAMPLE, ".bed")
            fwrite(bed_df, file = roh_bed_file, sep = "\t", quote = FALSE, row.names = FALSE)
          } # end sample
        } else {
          cat(sprintf("      File not found for Cycle %d, Run %d: %s\n", cycle, run, roh_file))
        } # end error from file check
      } # end cycle
    } # end run
    
    # Calculate mean and confidence intervals per chromosome, simulation, cycle, and sample
    #froh_summary <- froh_data %>%
    #  group_by(Chromosome, Simulation, Cycle, Sample) %>%
    #  summarise(
    #    Mean_FROH_2Mb = mean(FROH_2Mb, na.rm = TRUE),
    #    CI_Lower_2Mb = Mean_FROH_2Mb - qt(0.975, df = n() - 1) * sd(FROH_2Mb, na.rm = TRUE) / sqrt(n()),
    #    CI_Upper_2Mb = Mean_FROH_2Mb + qt(0.975, df = n() - 1) * sd(FROH_2Mb, na.rm = TRUE) / sqrt(n()),
    #    Mean_FROH_50Kb = mean(FROH_50Kb, na.rm = TRUE),
    #    CI_Lower_50Kb = Mean_FROH_50Kb - qt(0.975, df = n() - 1) * sd(FROH_50Kb, na.rm = TRUE) / sqrt(n()),
    #    CI_Upper_50Kb = Mean_FROH_50Kb + qt(0.975, df = n() - 1) * sd(FROH_50Kb, na.rm = TRUE) / sqrt(n()),
    #    Mean_FROH_all = mean(FROH_all, na.rm = TRUE),
    #    CI_Lower_all = Mean_FROH_all - qt(0.975, df = n() - 1) * sd(FROH_all, na.rm = TRUE) / sqrt(n()),
    #    CI_Upper_all = Mean_FROH_all + qt(0.975, df = n() - 1) * sd(FROH_all, na.rm = TRUE) / sqrt(n()),
    #    .groups = "drop" # Prevent unwanted grouping in the output
    #  )

    # Append to final results
    final_results <- rbind(final_results, froh_data)
  } # end simulation
} # end chromosome

final_results <- final_results %>%
  mutate(Simulation = factor(Simulation, levels = c("R", "N", "W", "D")))

# Write the final results to a file
write.table(final_results, file="FROH_Summary.txt", sep="\t", row.names=FALSE, quote=FALSE)

############################### PLOT

# Reshape data (wide → long)
library(dplyr)
library(tidyr)
library(ggplot2)

froh_long <- final_results %>%
  pivot_longer(
    cols = starts_with("FROH_"),
    names_to = "Scale",
    values_to = "FROH"
  )

# Compute per-run means
# collapses individuals → gives each run a single value per cycle.

run_means <- froh_long %>%
  group_by(Simulation, Cycle, Scale, Run) %>%
  summarise(
    Run_Mean = mean(FROH, na.rm = TRUE),
    .groups = "drop"
  )

#Compute between-run variance
between_run <- run_means %>%
  group_by(Simulation, Cycle, Scale) %>%
  summarise(
    Mean_FROH = mean(Run_Mean, na.rm = TRUE),
    SD_between = sd(Run_Mean, na.rm = TRUE),
    n_runs = n(),
    SE_between = SD_between / sqrt(n_runs),
    CI_Lower = Mean_FROH - qt(0.975, df = n_runs - 1) * SE_between,
    CI_Upper = Mean_FROH + qt(0.975, df = n_runs - 1) * SE_between,
    .groups = "drop"
  )

# Plot
ggplot() +
  
  # Each run = one line (NO cross-scale connections!)
  geom_line(
    data = run_means,
    aes(
      x = Cycle,
      y = Run_Mean,
      group = interaction(Run, Scale, Simulation),
      color = Scale
    ),
    alpha = 0.25
  ) +
  
  # Mean across runs
  geom_line(
    data = between_run,
    aes(
      x = Cycle,
      y = Mean_FROH,
      group = interaction(Scale, Simulation),
      color = Scale
    ),
    linewidth = 1.2
  ) +
  
  # Between-run CI
  geom_ribbon(
    data = between_run,
    aes(
      x = Cycle,
      ymin = CI_Lower,
      ymax = CI_Upper,
      group = interaction(Scale, Simulation),
      fill = Scale
    ),
    alpha = 0.2
  ) +
  
  facet_wrap(~ Simulation) +
  theme_minimal() +
  labs(
    title = "Between-Run Variation in FROH Over Time",
    x = "Cycle",
    y = "Mean FROH"
  )

ggsave("Plots/FROH_scale_Run_variation.pdf", width = 8, height = 6)
ggsave("Plots/FROH_scale_Run_variation.png", width = 8, height = 6)
