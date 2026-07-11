#########################################################
## Stacked

library(ggplot2)
library(dplyr)
library(readxl)
library(patchwork)
library(ggrepel)

# Load the data from Excel
#data <- read_excel("Ne_log_All.xlsx")
data <- read.csv("compiled_data.csv")
data <- data[, -ncol(data)]

colnames(data) <- c("Sim", "Run", "cycle", "N(t-1)", "N(t)", "Ne_heterozygosity", "heterozygosity")

# Compute summary statistics: mean and 95% CI for each simulation type per cycle
summary_data_Nt <- data %>%
  group_by(Sim, cycle) %>%
  summarise(
    mean_value = mean(`N(t)`, na.rm = TRUE),
    lower_CI = quantile(`N(t)`, 0.025, na.rm = TRUE),
    upper_CI = quantile(`N(t)`, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

summary_data_het <- data %>%
  group_by(Sim, cycle) %>%
  summarise(
    mean_value = mean(heterozygosity, na.rm = TRUE),
    lower_CI = quantile(heterozygosity, 0.025, na.rm = TRUE),
    upper_CI = quantile(heterozygosity, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot function for N(t) and heterozygosity with distinct y-axes
plot_sim_Nt <- function(sim_type) {
  ggplot() +
    # Plot N(t)
    #geom_line(data = filter(data, Sim == sim_type), 
    #          aes(x = cycle, y = `N(t)`, group = Run), 
    #          color = "blue", alpha = 0.1) +  
    # 95% CI shading for N(t)
    geom_ribbon(data = filter(summary_data_Nt, Sim == sim_type), 
                aes(x = cycle, ymin = lower_CI, ymax = upper_CI), 
                fill = "blue", alpha = 0.2) +
    # Mean line for N(t)
    geom_line(data = filter(summary_data_Nt, Sim == sim_type), 
              aes(x = cycle, y = mean_value), 
              color = "blue", linewidth = 1) +
    labs(x = "Cycle", y = "Population Size N(t)") +
    scale_y_continuous(limits = c(0, 510)) +
    theme_minimal() +
    theme(axis.title.y.left = element_text(color = "black"),
          axis.text.y.left = element_text(color = "black")) +
    # Add labels below the line at cycles 1, 10, 20
    geom_text_repel(data = filter(summary_data_Nt, Sim == sim_type, cycle %in% c(1, 10, 20)),
              aes(x = cycle, y = mean_value, label = round(mean_value, 0)),
              color = "black", size = 3, vjust = 1.5) +
    # Add points at cycles 1, 10, 20
    geom_point(data = filter(summary_data_Nt, Sim == sim_type, cycle %in% c(1, 10, 20)),
               aes(x = cycle, y = mean_value),
               color = "black", size = 1)
}

plot_sim_Het <- function(sim_type) {
  ggplot() +
    # Plot heterozygosity
    #geom_line(data = filter(data, Sim == sim_type), 
    #          aes(x = cycle, y = heterozygosity, group = Run), 
    #          color = "green", alpha = 0.1) +  
    # 95% CI shading for heterozygosity
    geom_ribbon(data = filter(summary_data_het, Sim == sim_type), 
                aes(x = cycle, ymin = lower_CI, ymax = upper_CI), 
                fill = "green", alpha = 0.2) +
    # Mean line for heterozygosity
    geom_line(data = filter(summary_data_het, Sim == sim_type), 
              aes(x = cycle, y = mean_value), 
              color = "green", size = 1) +
    labs(x = "Cycle", y = "Heterozygosity") +
    scale_y_continuous(limits = c(min(data$heterozygosity, na.rm = TRUE), 
                                  max(data$heterozygosity, na.rm = TRUE))) +
    theme_minimal() +
    theme(axis.title.y.left = element_text(color = "black"),
          axis.text.y.left = element_text(color = "black")) +
    # Add labels below the line at cycles 1, 10, 20
    geom_text_repel(data = filter(summary_data_het, Sim == sim_type, cycle %in% c(1, 10, 20)),
              aes(x = cycle, y = mean_value, label = round(mean_value, 7)),
              color = "black", size = 3, vjust = -1.5) +
    # Add points at cycles 1, 10, 20
    geom_point(data = filter(summary_data_het, Sim == sim_type, cycle %in% c(1, 10, 20)),
               aes(x = cycle, y = mean_value),
               color = "black", size = 1)
}

# Generate plots for each simulation type (N and W)
plot_Nt_N <- plot_sim_Nt("N")
plot_Het_N <- plot_sim_Het("N")
plot_Nt_W <- plot_sim_Nt("W")
plot_Het_W <- plot_sim_Het("W")
plot_Nt_R <- plot_sim_Nt("R")
plot_Het_R <- plot_sim_Het("R")

# Combine plots using patchwork
(plot_Nt_R | plot_Nt_N | plot_Nt_W) / (plot_Het_R | plot_Het_N | plot_Het_W)

# Save the combined plot as an SVG
ggsave("Plots/combined_3_stacked_plot.pdf", 
       plot = (plot_Nt_R | plot_Nt_N | plot_Nt_W) / (plot_Het_R | plot_Het_N | plot_Het_W), 
       device = "pdf", 
       width = 10, height = 5, dpi = 300)

###############################################################################
## Statistical Tests

library(lme4)
library(lmerTest)
library(emmeans)

nt_model <- lmer(
  `N(t)` ~ Sim * cycle + (1 | Run),
  data = data
)

cat("\n================ N(t) MODEL ================\n")
print(anova(nt_model))
print(summary(nt_model))

cat("\n---- N(t) slope comparisons ----\n")
print(emtrends(nt_model, pairwise ~ Sim, var = "cycle"))

het_model <- lmer(
  heterozygosity ~ Sim * cycle + (1 | Run),
  data = data
)

cat("\n================ HET MODEL (STACKED DATA) ================\n")
print(anova(het_model))
print(summary(het_model))

cat("\n---- HET slope comparisons ----\n")
print(emtrends(het_model, pairwise ~ Sim, var = "cycle"))


################################################################################
### Combined Line Plots with All Metrics

# Order simulations
summary_data_Nt$Sim  <- factor(summary_data_Nt$Sim, levels = c("R", "N", "W", "D"))
summary_data_het$Sim <- factor(summary_data_het$Sim, levels = c("R", "N", "W", "D"))

summary_data_Nt$cycle[summary_data_Nt$cycle == 1] <- 0
summary_data_het$cycle[summary_data_het$cycle == 1] <- 0

# Population size plot
p_Nt <- ggplot(
  summary_data_Nt,
  aes(x = cycle, y = mean_value, group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Cycle", y = "Population Size N(t)") +
  scale_y_continuous(limits = c(0, 510)) +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "black"),
    axis.text.y.left  = element_text(color = "black")
  ) +
  geom_point(
    data = filter(summary_data_Nt, cycle %in% c(1, 10, 20)),
    aes(x = cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_Nt, cycle %in% c(1, 10, 20)),
    aes(label = round(mean_value, 0)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) + theme(legend.position = "none")
p_Nt

# Heterozygosity plot
p_Het <- ggplot(
  summary_data_het,
  aes(x = cycle, y = mean_value, group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Cycle", y = "Heterozygosity") +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "black"),
    axis.text.y.left  = element_text(color = "black")
  ) +
  geom_point(
    data = filter(summary_data_het, cycle %in% c(1, 10, 20)),
    aes(x = cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_het, cycle %in% c(1, 10, 20)),
    aes(label = round(mean_value, 7)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) + theme_minimal() + theme(legend.position = "none") 
p_Het

### --- ADD IDrisk --- REQUIRES merged_results from 6.0.SLiM_plot_ID_risk--FINAL.R

library(dplyr)

summary_data_IDrisk <- merged_results %>%
  group_by(Simulation, Cycle) %>%
  summarise(
    mean_value = mean(ID_risk_2mb, na.rm = TRUE),
    sd_value   = sd(ID_risk_2mb, na.rm = TRUE),
    n          = n(),
    
    se         = sd_value / sqrt(n),
    
    lower_CI   = mean_value - 1.96 * se,
    upper_CI   = mean_value + 1.96 * se,
    
    .groups = "drop"
  )
colnames(summary_data_IDrisk)
unique(summary_data_IDrisk$Simulation)
unique(summary_data_IDrisk$Cycle)

summary_data_IDrisk$Simulation <- factor(summary_data_IDrisk$Simulation, levels = c("R", "N", "W", "D"))
summary_data_IDrisk$Sim <- summary_data_IDrisk$Simulation

p_IDrisk <-ggplot(
  summary_data_IDrisk,
  aes(x = Cycle, y = mean_value, group = Sim, color = Sim, fill = Sim)) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = filter(summary_data_IDrisk, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  labs(x = "Cycle", y = "mean IDrisk") +
  geom_text_repel(
    data = filter(summary_data_IDrisk, Cycle %in% c(0, 10, 20)),
    aes(label = round(mean_value, 3)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) + theme_minimal() + theme(legend.position = "none")
p_IDrisk

### --- ADD KING Kinship --- REQUIRES sample_summary from 7.1.plot_ngsRelate.R

summary_data_rel_king <- sample_summary %>%
  group_by(Simulation, Cycle) %>%
  #filter(Cycle %in% c(0,5,10,15,20)) %>%
  summarise(
    mean_value = mean(mean_KING, na.rm = TRUE),
    sd_value   = sd(mean_KING, na.rm = TRUE),
    n          = n(),
    
    se         = sd_value / sqrt(n),
    
    lower_CI   = mean_value - 1.96 * se,
    upper_CI   = mean_value + 1.96 * se,
    
    .groups = "drop"
  )

summary_data_rel_king$Sim <- factor(
  summary_data_rel_king$Simulation,
  levels = c("R", "N", "W", "D")
)

p_Rel_king <- ggplot(
  summary_data_rel_king,
  aes(x = Cycle, y = mean_value,
      group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = filter(summary_data_rel_king, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_rel_king, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value, label = sprintf("%.2f", mean_value)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(x = "Cycle", y = "Mean Kinship (KING)") + theme_minimal() + theme(legend.position = "none")
p_Rel_king

### --- ADD rab Kinship --- REQUIRES sample_summary from 7.1.plot_ngsRelate.R

summary_data_rel_rab <- sample_summary %>%
  group_by(Simulation, Cycle) %>%
  #filter(Cycle %in% c(0,5,10,15,20)) %>%
  summarise(
    mean_value = mean(mean_F, na.rm = TRUE),
    sd_value   = sd(mean_F, na.rm = TRUE),
    n          = n(),
    
    se         = sd_value / sqrt(n),
    
    lower_CI   = mean_value - 1.96 * se,
    upper_CI   = mean_value + 1.96 * se,
    
    .groups = "drop"
  )

summary_data_rel_rab$Sim <- factor(
  summary_data_rel_rab$Simulation,
  levels = c("R", "N", "W", "D")
)

p_Rel_rab <- ggplot(
  summary_data_rel_rab,
  aes(x = Cycle, y = mean_value,
      group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = filter(summary_data_rel_rab, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_rel_rab, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value, label = sprintf("%.2f", mean_value)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(x = "Cycle", y = "Mean Kinship (rab)") + theme_minimal() + theme(legend.position = "none")
p_Rel_rab

### --- ADD F1 Kinship --- REQUIRES sample_summary from 7.1.plot_ngsRelate.R

summary_data_rel_F1 <- merged_results %>%
  group_by(Simulation, Cycle) %>%
  #filter(Cycle %in% c(0,5,10,15,20)) %>%
  summarise(
    mean_value = mean(mean_F, na.rm = TRUE),
    sd_value   = sd(mean_F, na.rm = TRUE),
    n          = n(),
    
    se         = sd_value / sqrt(n),
    
    lower_CI   = mean_value - 1.96 * se,
    upper_CI   = mean_value + 1.96 * se,
    
    .groups = "drop"
  )

summary_data_rel_F1$Sim <- factor(
  summary_data_rel_F1$Simulation,
  levels = c("R", "N", "W", "D")
)

p_Rel_F1 <- ggplot(
  summary_data_rel_F1,
  aes(x = Cycle, y = mean_value,
      group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = filter(summary_data_rel_F1, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_rel_F1, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value, label = sprintf("%.2f", mean_value)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(x = "Cycle", y = "Mean Inbreeding (F)") + theme_minimal() + theme(legend.position = "none")
p_Rel_F1

### --- ADD theta Kinship --- REQUIRES sample_summary from 7.1.plot_ngsRelate.R

summary_data_rel_theta <- sample_summary %>%
  group_by(Simulation, Cycle) %>%
  #filter(Cycle %in% c(0,5,10,15,20)) %>%
  summarise(
    mean_value = mean(mean_theta, na.rm = TRUE),
    sd_value   = sd(mean_theta, na.rm = TRUE),
    n          = n(),
    
    se         = sd_value / sqrt(n),
    
    lower_CI   = mean_value - 1.96 * se,
    upper_CI   = mean_value + 1.96 * se,
    
    .groups = "drop"
  )

summary_data_rel_theta$Sim <- factor(
  summary_data_rel_theta$Simulation,
  levels = c("R", "N", "W", "D")
)

p_Rel_theta <- ggplot(
  summary_data_rel_theta,
  aes(x = Cycle, y = mean_value,
      group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = filter(summary_data_rel_theta, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_rel_theta, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value, label = sprintf("%.2f", mean_value)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(x = "Cycle", y = "Mean Kinship (theta)") + theme_minimal() + theme(legend.position = "none")
p_Rel_theta


### --- ADD FROH > 1Mb --- REQUIRES merged_results from 6.0.SLiM_plot_ID_risk--edit4FINAL.R

colnames(merged_results)

summary_data_FROH <- merged_results %>%
  group_by(Simulation, Cycle) %>%
  #filter(Cycle %in% c(0,5,10,15,20)) %>%
  summarise(
    mean_value = mean(FROH_1Mb, na.rm = TRUE),
    sd_value   = sd(FROH_1Mb, na.rm = TRUE),
    n          = n(),
    
    se         = sd_value / sqrt(n),
    
    lower_CI   = mean_value - 1.96 * se,
    upper_CI   = mean_value + 1.96 * se,
    
    .groups = "drop"
  )

summary_data_FROH$Sim <- factor(
  summary_data_FROH$Simulation,
  levels = c("R", "N", "W", "D")
)

p_FROH <- ggplot(
  summary_data_FROH,
  aes(x = Cycle, y = mean_value,
      group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = filter(summary_data_FROH, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_FROH, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value, label = sprintf("%.2f", mean_value)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(x = "Cycle", y = "Mean FRoH") + theme_minimal() + theme(legend.position = "none")
p_FROH

### --- ADD FROH > 1Mb --- REQUIRES merged_results from 6.0.SLiM_plot_ID_risk--edit4FINAL.R

colnames(merged_results)

summary_data_FROH2 <- merged_results %>%
  group_by(Simulation, Cycle) %>%
  #filter(Cycle %in% c(0,5,10,15,20)) %>%
  summarise(
    mean_value = mean(FROH_2Mb, na.rm = TRUE),
    sd_value   = sd(FROH_2Mb, na.rm = TRUE),
    n          = n(),
    
    se         = sd_value / sqrt(n),
    
    lower_CI   = mean_value - 1.96 * se,
    upper_CI   = mean_value + 1.96 * se,
    
    .groups = "drop"
  )

summary_data_FROH2$Sim <- factor(
  summary_data_FROH2$Simulation,
  levels = c("R", "N", "W", "D")
)

p_FROH2 <- ggplot(
  summary_data_FROH2,
  aes(x = Cycle, y = mean_value,
      group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = filter(summary_data_FROH2, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_FROH2, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value, label = sprintf("%.2f", mean_value)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(x = "Cycle", y = "Mean FRoH") + theme_minimal() + theme(legend.position = "none")
p_FROH2

### --- ADD HET per KB --- REQUIRES merged_results from 6.0.SLiM_plot_ID_risk--edit4FINAL.R

summary_data_HpK <- merged_results %>%
  group_by(Simulation, Cycle) %>%
  #filter(Cycle %in% c(0,5,10,15,20)) %>%
  summarise(
    mean_value = mean(HET_PER_KB, na.rm = TRUE),
    sd_value   = sd(HET_PER_KB, na.rm = TRUE),
    n          = n(),
    
    se         = sd_value / sqrt(n),
    
    lower_CI   = mean_value - 1.96 * se,
    upper_CI   = mean_value + 1.96 * se,
    
    .groups = "drop"
  )

summary_data_HpK$Sim <- factor(
  summary_data_HpK$Simulation,
  levels = c("R", "N", "W", "D")
)

p_HpK <- ggplot(
  summary_data_HpK,
  aes(x = Cycle, y = mean_value,
      group = Sim, color = Sim, fill = Sim)
) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = filter(summary_data_HpK, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value),
    inherit.aes = TRUE,
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = filter(summary_data_HpK, Cycle %in% c(0, 10, 20)),
    aes(x = Cycle, y = mean_value, label = sprintf("%.2f", mean_value)),
    color = "black",
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(x = "Cycle", y = "Mean HET / KB") + theme_minimal() + theme(legend.position = "none")
p_HpK

# List of Plots
p_Nt
p_IDrisk
p_Het
p_HpK
p_FROH
p_Rel_F1
p_Rel_theta
p_Rel_rab
p_Rel_king

# Combine plots
plot <- p_IDrisk / (p_Het + p_Rel_F1 + p_Rel_theta + p_Nt ) +
  plot_layout(heights = c(1, 2))
plot

# FINAL plot
plot <- p_IDrisk + p_Nt + p_FROH + p_Het + p_Rel_theta + p_Rel_king +
plot_layout(
  ncol = 2,
  nrow = 3,
  widths = 3,
  heights =2,
  guides = "collect")
plot

# Save the combined plot as an SVG
ggsave("Plots/Nt_Het_IDrisk_Rel.pdf", plot, width = 10, height = 10)
ggsave("Plots/Nt_Het_IDrisk_Rel.png", plot, width = 10, height = 10)
