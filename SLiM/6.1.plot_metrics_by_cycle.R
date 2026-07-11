library(ggplot2)
library(dplyr)
library(readr)
library(broom)

##################### Plot 1: IDrisk FRoH per Cycle

# Load FROH summary data
froh_all <- read_tsv("FROH_Summary.txt")
froh_zero <- read_tsv("ZERO/FROH_Summary.txt")

# Combine data
final_results <- bind_rows(froh_all, froh_zero)

# Summarize mean FROH per Sample-Cycle-Simulation
froh_summary_plot <- final_results %>%
  group_by(Sample, Cycle, Simulation, Run) %>%
  summarise(
    Mean = mean(FROH_2Mb, na.rm = TRUE),
  ) %>%
  ungroup()

# Focus on cycles of interest
target_cycles <- c(0, 5, 10, 15, 20)
froh_subset <- froh_summary_plot %>%
  filter(Cycle %in% target_cycles)
head(froh_subset)

# Set specific order
froh_subset$Simulation <- factor(froh_subset$Simulation, levels = c("R", "N", "W", "D"))

# --------------------------
# 1. Within-simulation tests: Does FROH change from 0 to 20?
# --------------------------
within_results <- froh_subset %>%
  filter(Cycle %in% c(0, 20)) %>%
  group_by(Simulation) %>%
  group_map(~ {
    mod <- lm(Mean ~ as.factor(Cycle), data = .x)
    broom::tidy(anova(mod))
  }, .keep = TRUE)

names(within_results) <- unique(froh_subset$Simulation)
print("ANOVA for FROH change from Cycle 0 to 20 (by simulation):")
print(within_results)

# --------------------------
# 2. Across-simulation test: Is the change in FROH different by simulation?
# --------------------------

# Filter to cycles 0 and 20
change_data <- froh_subset %>%
  filter(Cycle %in% c(0, 20))

# Fit model with interaction
interaction_model <- aov(Mean ~ Simulation * as.factor(Cycle), data = change_data)
summary(interaction_model)

# Optional: Tukey HSD post-hoc test
TukeyHSD(interaction_model)

# --------------------------
# 3. Visualization
# --------------------------
ggplot(froh_subset, aes(x = as.factor(Cycle), y = Mean)) +
  geom_violin(aes(fill = Simulation), alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  labs(title = "FROH (>2Mb) over Cycles by Simulation",
       x = "Cycle", y = "Mean FROH per individual") +
  facet_wrap(~ Simulation) +
  theme_minimal()

ggsave("FROH_over_cycles_by_sim.violin.pdf", height = 10, width = 30)

ggplot(froh_subset, aes(x = as.factor(Cycle), y = Mean, color = Simulation, group = Simulation)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(title = "Change in FROH over Time by Simulation",
       x = "Cycle",
       y = "Mean FROH (>2Mb)",
       color = "Simulation") +
  theme_minimal(base_size = 14)

ggsave("FROH_over_cycles_by_sim.stat_summary.pdf")

##################### Plot 2: IDrisk heterozygosity outide ROH per Cycle

# Load heterozygosity data
het <- read_csv("heterozygosity_per_kb_all_samples.csv")
het_zero <- read_csv("ZERO/heterozygosity_per_kb_all_samples.csv")

# Combine data
het_all <- bind_rows(het, het_zero)

# Summarize mean HET_PER_KB per Sample-Cycle-Simulation
het_summary <- het_all %>%
  group_by(Sample, Cycle, Simulation) %>%
  summarise(
    Mean = mean(HET_PER_KB, na.rm = TRUE)
  ) %>%
  ungroup()

# Focus on target cycles
target_cycles <- c(0, 5, 10, 15, 20)
het_subset <- het_summary %>%
  filter(Cycle %in% target_cycles)

# Set specific order
het_subset$Simulation <- factor(het_subset$Simulation, levels = c("R", "N", "W", "D"))

# --------------------------
# 1. Within-simulation tests: Does heterozygosity change from 0 to 20?
# --------------------------
within_het_results <- het_subset %>%
  filter(Cycle %in% c(0, 20)) %>%
  group_by(Simulation) %>%
  group_map(~ {
    mod <- lm(Mean ~ as.factor(Cycle), data = .x)
    broom::tidy(anova(mod))
  }, .keep = TRUE)

names(within_het_results) <- unique(het_subset$Simulation)
print("ANOVA for heterozygosity change from Cycle 0 to 20 (by simulation):")
print(within_het_results)

# --------------------------
# 2. Across-simulation test: Is the change in heterozygosity different by simulation?
# --------------------------

change_data_het <- het_subset %>%
  filter(Cycle %in% c(0, 20))

interaction_model_het <- aov(Mean ~ Simulation * as.factor(Cycle), data = change_data_het)
summary(interaction_model_het)

# Optional post-hoc
TukeyHSD(interaction_model_het)

# --------------------------
# 3. Visualization
# --------------------------

# Violin plot
ggplot(het_subset, aes(x = as.factor(Cycle), y = Mean)) +
  geom_violin(aes(fill = Simulation), alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  labs(title = "Heterozygosity per Kb over Cycles by Simulation",
       x = "Cycle", y = "Mean Heterozygosity (per Kb)") +
  facet_wrap(~ Simulation) +
  theme_minimal()

ggsave("HET_PER_KB_over_cycles_by_sim.violin.pdf", height = 10, width = 30)

# Line plot with error bars
ggplot(het_subset, aes(x = as.factor(Cycle), y = Mean, color = Simulation, group = Simulation)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(title = "Change in Heterozygosity per Kb over Time by Simulation",
       x = "Cycle",
       y = "Mean Heterozygosity (per Kb)",
       color = "Simulation") +
  theme_minimal(base_size = 14)

ggsave("HET_PER_KB_over_cycles_by_sim.stat_summary.pdf")



################################################################################
# ChatGPT aided version (tennis elbow day)

############################
# LOAD LIBRARIES
############################
library(ggplot2)
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(emmeans)

############################
# SETTINGS
############################
target_cycles <- c(0, 5, 10, 15, 20)

############################
# =========================
# 1. FROH ANALYSIS
# =========================
############################

# Load FROH summary data
froh_all <- read_tsv("FROH_Summary.txt")
froh_zero <- read_tsv("ZERO/FROH_Summary.txt")

# Combine data
final_results <- bind_rows(froh_all, froh_zero)

# Summarize mean FROH per Sample-Cycle-Simulation
froh_summary_plot <- final_results %>%
  group_by(Sample, Cycle, Simulation, Run) %>%
  summarise(
    Mean = mean(FROH_2Mb, na.rm = TRUE),
  ) %>%
  ungroup()

# Summarize FROH (KEEP Run + Sample)
froh_subset <- final_results %>%
  group_by(Sample, Cycle, Simulation, Run) %>%
  summarise(
    Mean = mean(FROH_2Mb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Cycle %in% target_cycles)

# Factor order
froh_subset$Simulation <- factor(froh_subset$Simulation, levels = c("R", "N", "W", "D"))

# --------------------------
# MIXED MODEL
# --------------------------
froh_model <- lmer(
  Mean ~ Simulation * Cycle + (1 | Run/Sample),
  data = froh_subset
)

cat("\n================ FROH MODEL ================\n")
print(anova(froh_model))
print(summary(froh_model))

# Pairwise slope differences (rate of change)
cat("\n---- FROH slope comparisons ----\n")
print(emtrends(froh_model, pairwise ~ Simulation, var = "Cycle"))

# --------------------------
# PLOTS
# --------------------------

# Violin + points
ggplot(froh_subset, aes(x = as.factor(Cycle), y = Mean)) +
  geom_jitter(aes(fill = Simulation, colour = Simulation), width = 0.1, alpha = 0.1, size = 1) +
  geom_violin(aes(fill = Simulation), alpha = 0.6, trim = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  facet_wrap(~ Simulation) +
  theme_minimal() +
  labs(title = "FROH (>2Mb) over Cycles by Simulation",
       x = "Cycle", y = "Mean FROH")

ggsave("Plots/FROH_violin.pdf", height = 10, width = 30)

froh_run <- froh_subset %>%
  group_by(Simulation, Run, Cycle) %>%
  summarise(
    Mean = mean(Mean, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(froh_run, aes(x = Cycle, y = Mean, color = Simulation)) +
  geom_line(aes(group = interaction(Simulation, Run)), alpha = 0.4, linewidth = 0.8) +
  stat_summary(
    aes(group = Simulation),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  stat_summary(
    aes(group = Simulation),
    fun = mean,
    geom = "point",
    size = 2
  ) +
  facet_wrap(~Simulation) +
  theme_minimal(base_size = 14) +
  labs(
    title = "FROH Trajectories (Run-level means)",
    x = "Cycle",
    y = "Mean FROH per Run"
  )

ggsave("Plots/FROH_trajectories.pdf")

############################
# =========================
# 2. HETEROZYGOSITY
# =========================
############################

het <- read_csv("heterozygosity_per_kb_all_samples.csv")
het_zero <- read_csv("ZERO/heterozygosity_per_kb_all_samples.csv")

het_all <- bind_rows(het, het_zero)

het_subset <- het_all %>%
  group_by(Sample, Cycle, Simulation, Run) %>%
  summarise(
    Mean = mean(HET_PER_KB, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Cycle %in% target_cycles)

het_subset$Simulation <- factor(het_subset$Simulation, levels = c("R", "N", "W", "D"))

# --------------------------
# MIXED MODEL
# --------------------------
het_model <- lmer(
  Mean ~ Simulation * Cycle + (1 | Run/Sample),
  data = het_subset
)

cat("\n================ HET MODEL ================\n")
print(anova(het_model))
print(summary(het_model))

cat("\n---- HET slope comparisons ----\n")
print(emtrends(het_model, pairwise ~ Simulation, var = "Cycle"))

# --------------------------
# PLOTS
# --------------------------

# Violin + points
ggplot(het_subset, aes(x = as.factor(Cycle), y = Mean)) +
  geom_jitter(aes(fill = Simulation, colour = Simulation), width = 0.1, alpha = 0.1, size = 1) +
  geom_violin(aes(fill = Simulation), alpha = 0.6, trim = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  facet_wrap(~ Simulation) +
  theme_minimal() +
  labs(title = "HETperKB over Cycles by Simulation",
       x = "Cycle", y = "Mean HET")

ggsave("Plots/HETperKB_violin.pdf", height = 10, width = 30)

het_run <- het_subset %>%
  group_by(Simulation, Run, Cycle) %>%
  summarise(
    Mean = mean(Mean, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(het_run, aes(x = Cycle, y = Mean, color = Simulation)) +
  geom_line(aes(group = interaction(Simulation, Run)), alpha = 0.4, linewidth = 0.8) +
  stat_summary(
    aes(group = Simulation),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  stat_summary(
    aes(group = Simulation),
    fun = mean,
    geom = "point",
    size = 2
  ) +
  facet_wrap(~Simulation) +
  theme_minimal(base_size = 14) +
  labs(
    title = "HET Trajectories (Run-level means)",
    x = "Cycle",
    y = "Mean HET per Run"
  )

ggsave("Plots/HETperKB_trajectories.pdf")

############################
# =========================
# 3. ID RISK
# =========================
############################

idrisk_subset <- merged_results %>%
  group_by(Sample, Cycle, Simulation, Run) %>%
  summarise(
    Mean = mean(ID_risk_2mb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Cycle %in% target_cycles)

idrisk_subset$Simulation <- factor(idrisk_subset$Simulation, levels = c("R", "N", "W", "D"))

# Violin plot
ggplot(idrisk_subset, aes(x = as.factor(Cycle), y = Mean)) +
  geom_violin(aes(fill = Simulation), alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  labs(title = "IDrisk over Cycles by Simulation",
       x = "Cycle", y = "Mean IDrisk") +
  facet_wrap(~ Simulation) +
  theme_minimal()

ggsave("Plots/IDrisk_over_cycles_by_sim.violin.pdf", height = 10, width = 30)

# Violin + points
ggplot(idrisk_subset, aes(x = as.factor(Cycle), y = Mean)) +
  geom_jitter(aes(fill = Simulation, colour = Simulation), width = 0.1, alpha = 0.1, size = 1) +
  geom_violin(aes(fill = Simulation), alpha = 0.6, trim = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  facet_wrap(~ Simulation) +
  theme_minimal() +
  labs(title = "IDrisk over Cycles by Simulation",
       x = "Cycle", y = "Mean IDrisk")

ggsave("Plots/IDrisk_violin.pdf", height = 10, width = 30)

idrisk_run <- idrisk_subset %>%
  group_by(Simulation, Run, Cycle) %>%
  summarise(
    Mean = mean(Mean, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(idrisk_run, aes(x = Cycle, y = Mean, color = Simulation)) +
  geom_line(aes(group = interaction(Simulation, Run)), alpha = 0.4, linewidth = 0.8) +
  stat_summary(
    aes(group = Simulation),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  stat_summary(
    aes(group = Simulation),
    fun = mean,
    geom = "point",
    size = 2
  ) +
  facet_wrap(~Simulation) +
  theme_minimal(base_size = 14) +
  labs(
    title = "IDrisk Trajectories (Run-level means)",
    x = "Cycle",
    y = "Mean IDrisk per Run"
  )

ggsave("Plots/IDrisk_trajectories.pdf")

############################
# 1. Prepare data (run-level consistency)
############################

idrisk_subset <- merged_results %>%
  group_by(Sample, Cycle, Simulation, Run) %>%
  summarise(
    Mean = mean(Mean_ID_risk_2mb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Cycle %in% target_cycles)

idrisk_subset$Simulation <- factor(idrisk_subset$Simulation,
                                   levels = c("R", "N", "W", "D"))

############################
# 2. YOUR ORIGINAL QUESTION: MEAN TRAJECTORY
############################

library(lme4)
library(lmerTest)

idrisk_model <- lmer(
  Mean ~ Simulation * Cycle + (1 | Run/Sample),
  data = idrisk_subset
)

cat("\n================ ID RISK MEAN MODEL ================\n")
print(anova(idrisk_model))
print(summary(idrisk_model))

############################
# 3. VISUAL: DISTRIBUTION (what you are noticing)
############################

ggplot(idrisk_subset, aes(x = Simulation, y = Mean, fill = Simulation)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3) +
  facet_wrap(~Cycle) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of Inbreeding Depression Risk",
       x = "Simulation",
       y = "Risk")

ggsave("IDrisk_distribution_by_simulation.pdf", height = 10, width = 14)

############################
# 4. QUANTIFY YOUR OBSERVATION: LOW-RISK PROPORTION
############################
# NOTE: adjust threshold if needed based on your scale

threshold <- quantile(idrisk_subset$Mean, 0.25, na.rm = TRUE)

low_risk_summary <- idrisk_subset %>%
  mutate(low_risk = Mean <= threshold) %>%
  group_by(Simulation, Cycle) %>%
  summarise(
    prop_low_risk = mean(low_risk),
    .groups = "drop"
  )

ggplot(low_risk_summary,
       aes(x = Cycle, y = prop_low_risk, color = Simulation)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  theme_minimal(base_size = 14) +
  labs(title = "Proportion of Low Inbreeding Depression Risk Individuals",
       y = "Proportion (lower quartile threshold)",
       x = "Cycle")

ggsave("Plots/IDrisk_low_risk_proportion.pdf", height = 6, width = 10)

############################
# 5. OPTIONAL: QUICK SUMMARY STATS
############################

idrisk_subset %>%
  group_by(Simulation, Cycle) %>%
  summarise(
    mean = mean(Mean),
    sd = sd(Mean),
    median = median(Mean),
    p10 = quantile(Mean, 0.1),
    p90 = quantile(Mean, 0.9),
    .groups = "drop"
  ) %>%
  print(n = 50)

