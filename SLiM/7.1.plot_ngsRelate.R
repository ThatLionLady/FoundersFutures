################################################################################
## ngsRelate of SLiM Output
################################################################################
library(tidyr)
library(ggplot2)

setwd("ZDC/SLiM/")

chr <- "A1"

# Define simulations 
simulations <- c("Reintro", "duoSupp", "noSupp", "w1Supp")
#simulations <- c("noSupp", "Reintro")

# Initialize list
rel_results_list <- list()

# To Test
#simulation <- "noSupp"
#run <- 0
#cycle <- 5

################################################################################
## COMPILE ngsRelate DATA w SLiM MATINGS
################################################################################
## PER SIMULATION 
################################################################################

for (simulation in simulations) {
  cat(sprintf("  Simulation: %s\n", simulation))
  
  sim_name <- toupper(substr(simulation, 1, 1))
  
  print(simulation)
  print(sim_name)
  
  ################################################################################
  ## PER RUN (gather matings.txt)
  ################################################################################
  
  for (run in 0:9) { # 10 runs labeled 0-9
    cat(sprintf("    Processing Run %d for Chromosome %s, Simulation %s\n", run, chr, sim_name))
    
    ################################################################################
    ## PER CYCLE
    ################################################################################
    
    for (cycle in 2:20) { # feature every 1 cycles for 20 cycles
      
      data <- paste0("A1_", sim_name, run, "_", cycle)
      file <- paste0("Relatedness/rel_", data, ".res")
      
      if (file.exists(file) && file.info(file)$size > 0) {
        cat(sprintf("    Processing Cycle %s\n", cycle))
        rel <- read.table(file ,h=T)
        # continue processing
      } else {
        cat(sprintf("    Skipping Cycle %s\n", cycle))
        next  # skip this iteration of the loop
      }
      
      # set IDs for merging
      # Make an ID map to rename default numbering
      id_map <- readLines(paste0("Lists/samples_", data, ".txt"))
      id_map
      
      # continue if TRUE
      max(rel$a, rel$b) + 1 <= length(id_map)
      
      rel$ida <- id_map[rel$a + 1]
      rel$idb <- id_map[rel$b + 1]
      unique(rel$ida)
      head(rel)
      
      rel$pair <- paste(
        pmin(rel$ida, rel$idb),
        pmax(rel$ida, rel$idb),
        sep = "-"
      )
      
      # Add SLiM run information
      rel$Simulation <- sim_name
      rel$Run <- run
      rel$Cycle <- cycle
      
      #######################
      # PEDIGREE INFERENCE
      #######################
      
      # Allele frequency-free inference of close familial relationships from genotypes or low depth sequencing data
      # https://doi.org/10.1101/260497
      ## "Since this table (the KING cutoffs) has overlapping kinship ranges for the PO and FS categories, 
      ## we used the R0 statistic to distinguish PO from FS relationships: Ignoring rare effects like germline
      ## mutations and genotyping errors the expected value for R0 for PO relatives is zero, while for FS the 
      ## value is above 0; we used an ad hoc cut-off of 0.02."
      # Following KING’s relationship inference framework (Manichaikul et al., 2010), 
      # first-degree relatives were classified using kinship coefficients, 
      # with parent–offspring and full-sibling relationships distinguished using the proportion of IBS0 sites (R0) (Waples et al., 2018)
      
      rel <- rel %>%
        mutate(
          relationship = case_when(
            KING > 0.354 ~ "Duplicate/MZ",
            
            # first-degree split
            #KING > 0.177 ~ "1st-degree",
            KING > 0.177 & R0 < 0.1 ~ "1st-dregree, Parent-Offspring (PO)",
            KING > 0.177 & R0 >= 0.1 ~ "1st-dregree, Full Siblings (FS)",
            
            # second-degree
            KING > 0.0884 ~ "2nd-degree",
            
            KING > 0.0442 ~ "3rd-degree",
            
            TRUE ~ "Unrelated"
          )
        )
      
      ## Append to results
      rel_results_list <- append(rel_results_list, list(rel))
      ################################################################################
    } #end cycle
  } #end run
}# end sim

# Combine all results into a single data.table
rel_combined <- rbindlist(rel_results_list)
colnames(rel_combined)

################################################################################
## ADD FOUNDERS (ZERO)
################################################################################

# Initialize list
zero_rel_results_list <- list()

for (simulation in simulations) {
  cat(sprintf("  Simulation: %s\n", simulation))
  
  sim_name <- toupper(substr(simulation, 1, 1))
  
  print(simulation)
  print(sim_name)
  
  ################################################################################
  ## PER RUN (gather matings.txt)
  ################################################################################
  
  cat(sprintf("    Processing Run %d for Chromosome %s, Simulation %s\n", run, chr, sim_name))
  
  ################################################################################
  ## PER CYCLE
  ################################################################################
  
  zero <- read.table(paste0("ZERO/rel_A1_", sim_name, "0_0.res") ,h=T)
  head(zero)
  
  # set IDs for merging
  # Make an ID map to rename default numbering
  id_map_list <- list(
    "19" = c(
      "KAR-F001","KAR-F002",
      "KHA-F001","KHA-F002",
      "MAK-F001","MAK-F002","MAK-F003","MAK-F004","MAK-F005","MAK-F006",
      "TEM-F001","TEM-F002","TEM-F003","TEM-F004",
      "TSW-F004","TSW-F005",
      "MOZ-M001",
      "TSW-M001","TSW-M002"
    ),
    "18" = c(
      "KAR-F001","KAR-F002",
      "KHA-F001","KHA-F002",
      "MAK-F001","MAK-F002","MAK-F003","MAK-F004","MAK-F005","MAK-F006",
      "TEM-F001","TEM-F002","TEM-F003","TEM-F004",
      "TSW-F004","TSW-F005",
      "TSW-M001","TSW-M002"
    )
  )
  
  id_map <- id_map_list[[as.character(length(unique(c(zero$a, zero$b))))]]
  id_map
  
  # continue if TRUE
  max(zero$a, zero$b) + 1 <= length(id_map)
  
  zero$ida <- id_map[zero$a + 1]
  zero$idb <- id_map[zero$b + 1]
  unique(zero$ida)
  head(zero)
  
  zero$pair <- paste(
    pmin(zero$ida, zero$idb),
    pmax(zero$ida, zero$idb),
    sep = "-"
  )
  
  # Add SLiM run information
  zero$Simulation <- sim_name
  zero$Run <- 0
  zero$Cycle <- 0
  
  #######################
  # PEDIGREE INFERENCE
  #######################
  
  zero <- zero %>%
    mutate(
      relationship = case_when(
        KING > 0.354 ~ "Duplicate/MZ",
        
        # first-degree split
        #KING > 0.177 ~ "1st-degree",
        KING > 0.177 & R0 < 0.1 ~ "1st-dregree, Parent-Offspring (PO)",
        KING > 0.177 & R0 >= 0.1 ~ "1st-dregree, Full Siblings (FS)",
        
        # second-degree
        KING > 0.0884 ~ "2nd-degree",
        
        KING > 0.0442 ~ "3rd-degree",
        
        TRUE ~ "Unrelated"
      )
    )
  
  zero_list <- lapply(0:9, function(run) {
    tmp <- zero
    tmp$Simulation <- sim_name
    tmp$Run <- run
    tmp$Cycle <- 0
    tmp
  })

  ## Append to results
  zero_rel_results_list <- append(zero_rel_results_list, zero_list)
  
}# end sim

# Combine all results into a single data.table
zero_combined <- rbindlist(zero_rel_results_list)
colnames(zero_combined)
head(zero_combined)
nrow(zero_combined)
zero_combined %>% distinct(across(all_of(c("Simulation", "Run", "Cycle"))))

rel_combined <- rbindlist(list(zero_combined, rel_combined), fill = TRUE)
head(rel_combined)

################################################################################
## CALCULATE MEAN PAIRWISE RELATEDNESS PER SAMPLE/SIM/RUN/CYCLE
################################################################################

# Make symmetric to include all combinations for each id as Sample
rel_sym <- rel_combined %>%
  transmute(
    Simulation, Run, Cycle,
    Sample = ida,
    Other  = idb,
    KING,
    theta,
    R0,
    R1,
    rab,
    FDiff,
    F1 = Fa,
    F2 = Fb
  ) %>%
  bind_rows(
    rel_combined %>%
      transmute(
        Simulation, Run, Cycle,
        Sample = idb,
        Other  = ida,
        KING,
        theta,
        R0,
        R1,
        rab,
        FDiff,
        F1 = Fb,
        F2 = Fa
      )
  )
head(rel_sym)

# get means across samples from the same simulation, run, cycle
sample_summary <- rel_sym %>%
  group_by(Simulation, Run, Cycle, Sample) %>%
  summarise(
    mean_KING = mean(KING, na.rm = TRUE),
    mean_theta = mean(theta, na.rm = TRUE),
    mean_rab = mean(rab, na.rm = TRUE),
    mean_F = mean(F1, na.rm = TRUE),
    mean_R0   = mean(R0, na.rm = TRUE),
    mean_R1   = mean(R1, na.rm = TRUE),
    n_pairs   = n(),
    .groups = "drop"
  )
head(sample_summary)
nrow(sample_summary %>% filter(Cycle == 0))
sample_summary %>% distinct(across(all_of(c("Simulation", "Run", "Cycle"))))

################################################################################
## COMBINE WITH FROH and HET
################################################################################

# Load FROH summary data
froh_all <- fread("FROH_Summary.txt")
froh_zero <- fread("ZERO/FROH_Summary.txt")
nrow(froh_zero)
froh_zero <- froh_zero %>%
  dplyr::select(-Run) %>%
  tidyr::crossing(Run = 0:9)
nrow(froh_zero)
froh_zero %>% distinct(across(all_of(c("Simulation", "Run", "Cycle"))))
final_results <- bind_rows(froh_all, froh_zero)
colnames(final_results)
nrow(final_results)

# Load HET summary data
het_all <- read.csv("heterozygosity_per_kb_all_samples.csv")
het_zero <- read.csv("ZERO/heterozygosity_per_kb_all_samples.csv")
het_zero <- het_zero %>%
  dplyr::select(-Run) %>%
  tidyr::crossing(Run = 0:9)
nrow(het_zero)
het_summary <- bind_rows(het_all, het_zero)
colnames(het_summary)
nrow(het_summary)

# Filter to check that data matches

sim_names <- toupper(substr(simulations, 1, 1))
keep_sims <- unique(sim_names)
keep_sims

sample_cycle_focus <- sample_summary %>% filter(Cycle %in% c("0", "5", "10", "15", "20"))
colnames(sample_cycle_focus)
nrow(sample_cycle_focus)
sample_cycle_focus %>% filter(Run == 0, Cycle == 0) %>% nrow()

filtered_froh <- final_results %>% filter(Simulation %in% keep_sims)
nrow(filtered_froh)
filtered_froh %>% filter(Run == 0, Cycle == 0) %>% nrow()

filtered_het <- het_summary %>% filter(Simulation %in% keep_sims)
nrow(filtered_het)
filtered_het %>% filter(Run == 0, Cycle == 0) %>% nrow()

keys <- c("Sample", "Simulation", "Run", "Cycle")

# Check that samples match across analyses:
count_results <- list()

for (sim_name in sim_names) {
  for (run in 0) { 
    for (cycle in c(0,5,10,15,20)) {
      het_count <- filtered_het %>% filter(Simulation == sim_name, Run == run, Cycle == cycle) %>% nrow()
      froh_count <- filtered_froh %>% filter(Simulation == sim_name, Run == run, Cycle == cycle) %>% nrow()
      rel_count <- sample_cycle_focus %>% filter(Simulation == sim_name, Run == run, Cycle == cycle) %>% nrow()
      
      count_results[[length(count_results) + 1]] <- data.frame(
        Simulation = sim_name,
        Run = run,
        Cycle = cycle,
        HET = het_count,
        FROH = froh_count,
        REL = rel_count)
    }
  }
}

count_results_df <- dplyr::bind_rows(count_results)
count_results_df

# IF THEY DON'T MATCH
# Get list of sample names and determine which ones don't match between 
sim_name = "R"
Run = 0
cycle = 5

rel_samples <- sample_summary %>%
  filter(Simulation == sim_name, Run == run, Cycle == cycle) %>%
  distinct(Sample) %>%
  arrange(Sample)

froh_samples <- final_results %>%
  filter(Simulation == sim_name, Run == run, Cycle == cycle) %>%
  distinct(Sample) %>%
  arrange(Sample)

het_samples <- het_summary %>%
  filter(Simulation == sim_name, Run == run, Cycle == cycle) %>%
  distinct(Sample) %>%
  arrange(Sample)

anti_join(rel_samples, froh_samples, by = "Sample")
anti_join(rel_samples, het_samples, by = "Sample")
anti_join(het_samples, froh_samples, by = "Sample")

# Merge FROH and heterozygosity results
# if there's an error for unused arguments
#detach("package:MASS", unload = TRUE)
library(dplyr)

merged_results <- filtered_froh %>%
  select(Sample, Simulation, Run, Cycle, FROH_2Mb, FROH_1Mb) %>%
  left_join(filtered_het %>%
              select(Sample, Simulation, Run, Cycle, HET_PER_KB),
            by = keys) %>%
  left_join(sample_cycle_focus %>%
              select(Sample, Simulation, Run, Cycle,
                     mean_theta, mean_KING, mean_rab, mean_F),
            by = keys)

# Calculate ID risk
merged_results <- merged_results %>%
  mutate(ID_risk_2mb = HET_PER_KB * FROH_2Mb)
colnames(merged_results)

# Ensure ID_risk_2mb is numeric
merged_results$ID_risk_2mb <- as.numeric(merged_results$ID_risk_2mb)

# Set levels
merged_results$Simulation <- factor(merged_results$Simulation, levels = c("R", "N", "W", "D"))
rel_combined$Simulation <- factor(rel_combined$Simulation, levels = c("R", "N", "W", "D"))

# Check that data is complete
merged_results[!complete.cases(merged_results), ] 

################################################################################
## Plot average pairwise relatedness by cycle

# per pair (mean pairwise relatedness per run across samples)
rel_combined %>%
  group_by(Simulation, Cycle, Run) %>%
  summarise(mean_KING = mean(KING, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(Cycle, mean_KING, color = Simulation)) +

  # individual runs
  geom_line(aes(group = interaction(Simulation, Run)), alpha = 0.3) +
  
  # simulation mean
  stat_summary(
    aes(group = Simulation),
    fun = mean,
    geom = "line",
    linewidth = 1.3
  )

ggsave("Plots/meanKINGxCycle.pdf", width = 10, height = 8)
ggsave("Plots/meanKINGxCycle.png", width = 10, height = 8)

# sample summary box plot (mean per sample pairwise relatedness per run) <-- results in the same plot! haha
sample_summary %>%
  group_by(Simulation, Run, Cycle) %>%
  summarise(mean_KING = mean(mean_KING, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(
    x = Cycle,
    y = mean_KING,
    color = Simulation,
    group = interaction(Simulation, Run)
  )) +
  geom_line(alpha = 0.5, linewidth = 1)

# Plot all metrics
library(tidyr)

plot_df <- merged_results %>%
  pivot_longer(
    cols = c(FROH_1Mb, HET_PER_KB, ID_risk_2mb, mean_KING,
             mean_rab, mean_F),
    names_to = "metric",
    values_to = "value"
  )

plot_df$metric <- factor(plot_df$metric,
                         levels = c("FROH_1Mb", "HET_PER_KB", "ID_risk_2mb",
                                    "mean_KING", "mean_rab", "mean_F"))

plot_df %>%
  group_by(Simulation, Run, Cycle, metric) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Cycle, y = value, color = Simulation)) +
  geom_line(aes(group = interaction(Simulation, Run)), alpha = 0.25) +
  stat_summary(aes(group = Simulation),
               fun = mean,
               geom = "line",
               linewidth = 1.2) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw()

ggsave("Plots/popgen_metrics.pdf", width = 10, height = 5)
ggsave("Plots/popgen_metrics.png", width = 10, height = 5)

################################################################################
## Relationships plot

rel_combined %>% 
  filter(KING > 0) %>% 
  filter(Run == 0) %>% 
  ggplot(aes(
    x = FDiff,
    y = KING,
    color = relationship
  )) +
  geom_point(alpha = 0.4, size = 1.5) +
  facet_grid(Simulation ~ Cycle) +
  theme_bw()

ggsave("Plots/KINGgt0_FDiff.pdf", width = 20, height = length(simulations))
ggsave("Plots/KINGgt0_FDiff.png", width = 20, height = length(simulations))

# King x R1 wrapped to see if the trend holds across cycles
rel_combined %>% 
  filter(KING > 0) %>%
  filter(Run == 0, Simulation == "N") %>%
  ggplot(aes(x = R1, y = KING, color = relationship)) +
  geom_point(alpha = 0.3, size = 1) +
  facet_wrap(~ Cycle) +
  theme_bw()

### HIGH NUMBER OF MZ/Dups?
# True inbreeding / drift in small populations (biological explanation)
# Pedigree compression effect (VERY important in simulations) > extreme IBD inflation due to pedigree collapse
# Threshold artefact > 0.354 is too low for this population because
#   *Inbred full siblings can sometimes exceed this
#   *highly bottlenecked populations can push many pairs above threshold
#   *genotyping error / low coverage can inflate KING upward
#   *Increasing counts can reflect “distribution shift into extreme relatedness tail” not literal duplication events

# Investigate! Plot all pairwise. Do the patterns between simulations and cycles change?
rel_combined %>% 
  filter(KING > 0) %>% 
  #filter(Run == 0) %>%
  ggplot(aes(Cycle, KING, color = Simulation)) +
  geom_jitter(alpha = 0.05, size = 0.5) +
  facet_wrap(~ Simulation, ncol = 1) 

ggsave("Plots/KINGxCycle_ALL.pdf", width = 20, height = length(simulations)*2)
ggsave("Plots/KINGxCycle_ALL.png", width = 20, height = length(simulations)*2)

# Look at proportion of high relatedness over time
rel_combined %>% 
  filter(Run == 0) %>% 
  group_by(Simulation, Cycle) %>% 
  summarise(
    prop_dup = mean(KING > 0.354, na.rm = TRUE),
    prop_1st = mean(KING > 0.177, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(prop_dup, prop_1st),
    names_to = "class",
    values_to = "prop"
  ) %>%
  ggplot(aes(Cycle, prop, color = class)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Simulation, ncol = 1) +
  theme_bw()

ggsave("Plots/proportion_high_rel.pdf", width = 10, height = length(simulations)*2)
ggsave("Plots/proportion_high_rel.png", width = 10, height = length(simulations)*2)

# Table of above
rel_combined %>%
  #filter(Run == 0) %>%
  group_by(Simulation, Cycle) %>%
  summarise(
    prop_dup = mean(KING > 0.354, na.rm = TRUE),
    prop_1st = mean(KING > 0.177, na.rm = TRUE),
    .groups = "drop"
  ) %>% print(n=80)

# RESULT for N:
#   prop_dup steadily declines: 0 → ~0.00075
#   prop_1st rises early, peaks, then declines: ~0 → 0.26 → 0.09
# early accumulation of relatedness → then dilution or restructuring of variation
# >> transient inbreeding signal followed by recovery / redistribution of diversity

# RESULT for R:
#   prop_dup: stays ~0 almost entirely, 0 → ~0.01 → ~0.001
# no strong MZ explosion
#   prop_1st: increases strongly: 0.01 → 0.28 (cycle 10) → 0.096
# ✔ clear peak in relatedness structure
# >> population becomes highly interconnected genetically

# Plot how the composition of relatedness changes over time
rel_classes <- rel_combined %>%
  filter(Run == 0) %>%
  mutate(
    class = case_when(
      KING > 0.354 ~ "Dup/MZ",
      KING > 0.177 ~ "1st",
      KING > 0.0884 ~ "2nd",
      KING > 0.0442 ~ "3rd",
      TRUE ~ "Unrelated"
    )
  )

rel_props <- rel_classes %>%
  group_by(Simulation, Cycle, class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Simulation, Cycle) %>%
  mutate(prop = n / sum(n))

rel_props$class <- factor(
  rel_props$class,
  levels = c("Dup/MZ", "1st", "2nd", "3rd", "Unrelated")
)

ggplot(rel_props, aes(Cycle, prop, fill = class)) +
  geom_area(position = "fill") + 
  facet_wrap(~ Simulation, ncol = 1) +
  theme_bw() 

ggsave("Plots/proportions_of_relationship_classes_over_time.pdf", width = 10, height = length(simulations)*2)
ggsave("Plots/proportions_of_relationship_classes_over_time.png", width = 10, height = length(simulations)*2)
