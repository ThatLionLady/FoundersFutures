library(ggplot2)
library(dplyr)
library(readr)
library(car)

# Featured cycle
cycle <- "20"

# Load the data
df <- read_tsv("counts4IDrisk.txt")

# Reorder Simulation factor
df$Simulation <- factor(df$Simulation, levels = c("R", "N", "W"))

# Make sure columns are in correct format
df <- df %>%
  mutate(
    Cycle = as.integer(Cycle),
    Run = as.factor(Run), # Treat run as categorical
    Chromosome = factor(Chromosome, levels = unique(Chromosome)) # Optional: maintain order
  )

# Levene’s test for each simulation (Testing if the variance in N_samples differs across chromosomes)
df %>% 
  filter(Cycle == cycle) %>% 
  group_by(Simulation) %>% 
  group_map(~ leveneTest(N_samples ~ Chromosome, data = .x))

# Mean differences (if certain chromosomes systematically have more or fewer samples than others within a simulation)
df %>%
  filter(Cycle == cycle) %>%
  group_by(Simulation) %>%
  group_map(~ aov(N_samples ~ Chromosome, data = .x) %>% summary())

# Plot
ggplot(df %>% filter(Cycle == cycle),
       aes(x = Chromosome, y = N_samples, fill = Simulation, shape = Simulation)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  theme_minimal(base_size = 10) +
  labs(title = paste0("Variance in Sample Counts by Chromosome (Cycle ",cycle,")"),
       y = "Number of Samples")

df %>%
  filter(Cycle == cycle) %>%
  ggplot(aes(x = Chromosome, y = N_samples, fill = Simulation)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ Simulation) +
  theme_minimal(base_size = 10) +
  labs(title = paste0("Variance in Sample Counts by Chromosome (Cycle ",cycle,")"),
       y = "Number of Samples")

################## Iterate!

cycles <- c(5,10,15,20)
df <- read_tsv("counts4IDrisk.txt") # Load the data
df$Simulation <- factor(df$Simulation, levels = c("R", "N", "W")) # Reorder Simulation factor

results <- data.frame(Test=character(), Cycle=integer(), Simulation=numeric(), `F-value`=numeric(), `p-value`=numeric(), Interpretation=character(), stringsAsFactors=FALSE)

for (cycle in cycles) {

  # Filter data for the selected cycle
  cycle_df <- df %>% filter(Cycle == cycle)
  
  # Run tests and collect results
  cycle_results <- cycle_df %>%
    group_by(Simulation) %>%
    group_map(~ {
      # Levene's test
      lev <- car::leveneTest(N_samples ~ factor(Chromosome), data = .x, center = "median")
      lev_f <- lev[1, "F value"]
      lev_p <- lev[1, "Pr(>F)"]
      lev_interp <- ifelse(lev_p < 0.05, "Significant variance", "Not significant")
      
      # ANOVA
      aov_result <- summary(aov(N_samples ~ Chromosome, data = .x))[[1]]
      aov_f <- aov_result[1, "F value"]
      aov_p <- aov_result[1, "Pr(>F)"]
      aov_interp <- ifelse(aov_p < 0.05, "Significant mean diff", "No significant diff")
      
      # Return one row of results
      tibble(
        Test = c("Levene's Test", "ANOVA"),
        Cycle = cycle,
        Simulation = unique(.y$Simulation),
        `F-value` = c(lev_f, aov_f),
        `p-value` = c(lev_p, aov_p),
        Interpretation = c(lev_interp, aov_interp)
      )
    }) %>%
    bind_rows()
  
  # Plot
  ggplot(cycle_df %>% filter(Cycle == cycle),
         aes(x = Chromosome, y = N_samples, fill = Simulation, shape = Simulation)) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    theme_minimal(base_size = 14) +
    labs(title = paste0("Variance in Sample Counts by Chromosome (Cycle ",cycle,")"),
         y = "Number of Samples")
  
  ggsave(paste0("variance_plot_cycle",cycle,".pdf"), width = 12, height = 6)
  
  cycle_df %>%
    filter(Cycle == cycle) %>%
    ggplot(aes(x = Chromosome, y = N_samples, fill = Simulation)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~ Simulation) +
    theme_minimal(base_size = 14) +
    labs(title = paste0("Variance in Sample Counts by Chromosome (Cycle ",cycle,")"),
         y = "Number of Samples")

  ggsave(paste0("variance_plot_by_model_cycle",cycle,".pdf"), width = 12, height = 6)
  
  results <- rbind(results, cycle_results)

}

# Write to file (TSV)
write.table(results, paste0("Plots/byCHR_statistical_results.tsv"), sep = "\t", row.names = FALSE)

