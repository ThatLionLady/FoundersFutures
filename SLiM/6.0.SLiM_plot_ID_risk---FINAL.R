library(data.table)
library(ggplot2)
library(MASS)
library(dplyr)

setwd("ZDC/SLiM")

# Load Inbreeding Depression Risk Data
idrisk_data <- read.csv("ROH_het_by_species.csv")
names(idrisk_data) <- c("infoID", "HET_PER_KB", "Mean_FROH_2Mb", "ID_risk_2mb")

# Load FROH summary data
froh_all <- fread("FROH_Summary.txt")
froh_zero <- fread("ZERO/FROH_Summary.txt")
final_results <- bind_rows(froh_all, froh_zero)
colnames(final_results)
# Load HET summary data
het_all <- read.csv("heterozygosity_per_kb_all_samples.csv")
het_zero <- read.csv("ZERO/heterozygosity_per_kb_all_samples.csv")
het_summary <- bind_rows(het_all, het_zero)
colnames(het_summary)

# Merge FROH and heterozygosity results
merged_results <- merge(final_results[, c("Sample", "FROH_2Mb", "Chromosome", "Simulation", "Cycle", "Run")],
  het_summary[, c("Sample", "Chromosome", "Simulation", "Cycle", "Run", "HET_PER_KB")],
  by = c("Sample", "Chromosome", "Simulation", "Cycle", "Run")
)

# Add infoID combining all information
merged_results$infoID <- paste(merged_results$Sample, merged_results$Simulation, 
                               merged_results$Cycle, merged_results$Run, sep = "_")

# Calculate ID risk
merged_results <- merged_results %>%
  mutate(ID_risk_2mb = HET_PER_KB * FROH_2Mb)
colnames(merged_results)

# Ensure ID_risk_2mb is numeric
merged_results$ID_risk_2mb <- as.numeric(merged_results$ID_risk_2mb)

################################## GENOMIC MEAN (across chromosomes) TO USE AS FINAL!!!!!!

######################################################
# STEP 1: Calculate Mean Across Chromosomes
######################################################

# Calculate mean values across chromosomes and runs for each Sample, Simulation, and Cycle
mean_across_chromosomes <- merged_results %>%
  group_by(Sample, Simulation, Cycle, Run) %>%
  summarise(
    Mean_ID_risk_2mb = mean(ID_risk_2mb, na.rm = TRUE),
    chrMean_FROH_2Mb = mean(FROH_2Mb, na.rm = TRUE),
    Mean_HET = mean(HET_PER_KB, na.rm = TRUE),
    n_runs = n_distinct(Run),
    n_chr = n_distinct(Chromosome),
    .groups = "drop"
  )

# Write the final results to a file
write.table(final_results, file="IDrisk_mean_across_chromosomes.txt", sep="\t", row.names=FALSE, quote=FALSE)

######################################################
# STEP 2: ggPlot Mean Values by Simulation and Cycle
######################################################

# Define a color palette for the plots
color_palette <- colorRampPalette(c("blue", "green", "orange", "red"))

# Create the plot
plot <- ggplot(mean_across_chromosomes, aes(x = chrMean_FROH_2Mb, y = Mean_HET, color = Mean_ID_risk_2mb)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradientn(colors = color_palette(100), name = "ID Risk (Mean)") +
  facet_grid(Simulation ~ Cycle) +
  theme_minimal() +
  labs(
    title = "Mean ID Risk Across Chromosomes by Simulation and Cycle",
    x = "Mean FROH (2Mb)",
    y = "Mean HET (non-ROH)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )
plot

# Save the plot to a file
ggsave("mean_id_risk_across_chromosomes.pdf", plot = plot, width = 12, height = 8)
ggsave("mean_id_risk_across_chromosomes.png", plot = plot, width = 12, height = 8)

######################################################
# STEP 2.a: Plot the Mean Data by Simulation and Cycle (baseR)
######################################################

# Create a data frame for risk lines
risk_lines <- data.frame(
  x = seq(0.01, 0.7, length.out = 100),
  vhigh_risk = 0.5 / seq(0.01, 0.7, length.out = 100),
  high_risk = 0.25 / seq(0.01, 0.7, length.out = 100),
  moderate_risk = 0.05 / seq(0.01, 0.7, length.out = 100)
)

# Define a color palette for ID risk
ryPal <- colorRampPalette(c('gold1', 'darkred'))

# Create a gradient color column for ID risk
# Define breaks that match your curves
breaks <- c(0, 0.05, 0.25, 0.5, Inf)

# Assign colors based on thresholds
mean_across_chromosomes$Col <- ryPal(length(breaks) - 1)[
  as.numeric(cut(mean_across_chromosomes$Mean_ID_risk_2mb,
                 breaks = breaks,
                 include.lowest = TRUE))
]
# Get unique combinations of Chromosome, Simulation, and Cycle
unique_combinations <- mean_across_chromosomes %>%
  dplyr::select(Simulation, Cycle) %>%
  dplyr::distinct()

# Loop through each combination and create a plot
for (i in seq_len(nrow(unique_combinations))) {
  # Subset the data for the current combination
  current_combination <- unique_combinations[i, ]
  sim <- current_combination$Simulation
  cycle <- current_combination$Cycle
  
  subset_data <- mean_across_chromosomes %>%
    filter(Simulation == sim, Cycle == cycle)
  
  focus_run <- subset_data %>% filter(Run == 0)
  other_runs <- subset_data %>% filter(Run != 0)
  
  # Define the output file name
  output_file <- paste0("Plots/Mean_IDrisk_2Mb_", sim, "_Cycle", cycle, ".pdf")
  
  # Create the plot
  pdf(output_file, width = 10, height = 8)
  
  # Scatter plot
  plot(other_runs$chrMean_FROH_2Mb,
       other_runs$Mean_HET,
       pch = 21,
       bg  = adjustcolor(other_runs$Col, alpha.f = 0.25),  # 25% opacity
       col = NA,
       lwd = 0.25,
       cex = 1,
       xlim = c(0, 0.6), ylim = c(0, 2.7),
       ylab = "Het/kb in non-ROH regions",
       xlab = "FROH",
       cex.axis = 1.1, cex.lab = 1.3,
       main = paste("Simulation:", sim, "Cycle:", cycle))
  
  # overlay focus run
  points(focus_run$chrMean_FROH_2Mb,
         focus_run$Mean_HET,
         pch = 21,
         bg  = focus_run$Col,
         col = "black",
         lwd = 0.8,
         cex = 1.5)
  
  # Add sample labels
  #text(subset_data$chrMean_FROH_2Mb, subset_data$Mean_HET,
  #     labels = subset_data$Sample, cex = 0.5, adj = c(0, 0.5), srt = -20)
  
  # Add risk lines
  x <- seq(0, 0.7, length.out = 100)
  vhigh_risk <- 0.5 / x
  high_risk <- 0.25 / x
  moderate_risk <- 0.05 / x
  lines(x, vhigh_risk, col = "darkred", lty = 2, lwd = 2)
  lines(x, high_risk, col = "red", lty = 2, lwd = 2)
  lines(x, moderate_risk, col = "orange", lty = 2, lwd = 2)
  
  # Add a legend to describe risk lines
  legend("topright", legend = c("Very High Risk", "High Risk", "Moderate Risk"),
         col = c("darkred", "red", "orange"), lty = 2, lwd = 2, bty = "n")
  
  # Add the color legend for ID risk
  colors <- ryPal(20)
  par(mar = c(5, 4, 4, 6)) # Adjust margins to make space for the legend
  image(1, seq(0, 0.9, length.out = 100), t(seq(0, 0.9, length.out = 100)),
        col = colors, axes = FALSE, xlab = "", ylab = "")
  axis(4, at = seq(0, 0.9, by = 0.1), labels = seq(0, 0.9, by = 0.1), las = 1)
  mtext("IDrisk = FROH * non-ROH Heterozygosity", side = 4, line = 2.5)
  
  # Close the PDF
  dev.off()
  
  ######################################################
  # STEP 2.b: Plot the Mean Data by Simulation and Cycle (ggplot)
  ######################################################
  
  # Plotting with ggplot2
  ggplot() +
    
    # background runs (faded)
    geom_point(
      data = other_runs,
      aes(x = chrMean_FROH_2Mb,
          y = Mean_HET,
          color = Mean_ID_risk_2mb),
      size = 3,
      alpha = 0.15
    ) +
    
    # focus run (strong)
    geom_point(
      data = focus_run,
      aes(x = chrMean_FROH_2Mb,
          y = Mean_HET,
          color = Mean_ID_risk_2mb),
      size = 3.5,
      alpha = 1
    ) +
    
    # risk gradient
    scale_color_gradient(low = "gold1", high = "darkred") +
    
    # risk threshold lines
    geom_line(data = risk_lines, aes(x = x, y = vhigh_risk), color = "red") +
    geom_line(data = risk_lines, aes(x = x, y = high_risk), color = "orange") +
    geom_line(data = risk_lines, aes(x = x, y = moderate_risk), color = "yellow") +
    
    # labels + theme
    labs(
      x = expression(F[ROH]),
      y = "Het/kb in non-ROH regions",
      color = "ID risk"
    ) +
    
    theme_minimal() +
    theme(
      legend.key.height = unit(1, "null"),
      legend.margin = margin(0, 0, 0, 0)
    ) +
    
    xlim(0, 0.6) +
    ylim(0, 2.7)
    
  
  ggsave(paste0("Plots/ggplot_Mean_IDrisk_2Mb_", sim, "_Cycle", cycle, ".pdf"))
  ggsave(paste0("Plots/ggplot_Mean_IDrisk_2Mb_", sim, "_Cycle", cycle, ".png"))
}

################################################################################

# Ensure factors are in order
mean_across_chromosomes <- mean_across_chromosomes %>%
  mutate(
    Simulation = factor(Simulation, levels = c("R", "N", "W", "D")),
    Cycle = factor(Cycle)
  )

# define focus vs background
focus_runs <- mean_across_chromosomes %>% filter(Run == 0)
other_runs <- mean_across_chromosomes %>% filter(Run != 0)

# Plot cycles across, sims down
ggplot() +
  
  # risk lines (same in every panel)
  geom_line(data = risk_lines, aes(x = x, y = vhigh_risk), color = "grey") +
  geom_line(data = risk_lines, aes(x = x, y = high_risk), color = "grey") +
  geom_line(data = risk_lines, aes(x = x, y = moderate_risk), color = "grey") +
    
  # run
  geom_point(
    data = mean_across_chromosomes,
    aes(x = chrMean_FROH_2Mb,
        y = Mean_HET,
        color = Mean_ID_risk_2mb),
    shape = 16,
    size = 2,
    alpha = 0.1
  ) +
  
  # risk gradient
  scale_color_gradient(low = "gold1", high = "darkred") +
  
  stat_density_2d(
    data = other_runs,
    aes(x = chrMean_FROH_2Mb,
        y = Mean_HET),
    size =3) +
  
  # faceting layout
  facet_grid(Simulation ~ Cycle) +
  
  labs(
    x = expression(F[ROH]),
    y = "Het/kb in non-ROH regions",
    color = "ID risk"
  ) +
  
  theme_minimal() +
  
  theme(
    strip.text = element_text(size = 10),
    legend.key.height = unit(1, "null"),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  
  xlim(0, 0.6) +
  ylim(0, 1.0)

ggsave("Plots/ggplot_IDrisk_Density.pdf", width = 10, height = 4)
ggsave("Plots/ggplot_IDrisk_Density.png", width = 10, height = 4)

###############################
## Density Plot (base R)

# turns the FROH x heterozygosity point cloud unto a continuous surface of point density, 
# contour lines represent elevation lines where mountains = high density of points and valleys = low density

make_plot <- function(bg_run, focus_run, sim, cycle) {
  
  # -------------------------
  # Base scatter
  # -------------------------
  plot(bg_run$chrMean_FROH_2Mb,
       bg_run$Mean_HET,
       pch = 21,
       bg  = adjustcolor(bg_run$Col, alpha.f = 0.15),
       col = NA,
       cex = 1 * facet_scale,
       xlim = c(0, 0.6),
       ylim = c(0, 2.7),
       xlab = expression(F[ROH]),
       ylab = "Het/kb in non-ROH regions",
       main = paste("Simulation:", sim, "Cycle:", cycle))
  
  
  # -------------------------
  # Focus overlay
  # -------------------------
  if (nrow(focus_run) > 0) {
    points(focus_run$chrMean_FROH_2Mb,
           focus_run$Mean_HET,
           pch = 21,
           bg  = focus_run$Col,
           col = "black",
           lwd = 1.2,
           cex = 1.4 * facet_scale)
  }
  
  # -------------------------
  # Density (independent of focus_run!)
  # -------------------------
  if (nrow(subset_run) >= 5 &&
      all(is.finite(subset_run$chrMean_FROH_2Mb)) &&
      all(is.finite(subset_run$Mean_HET))) {
    
    dens <- MASS::kde2d(subset_run$chrMean_FROH_2Mb,
                        subset_run$Mean_HET,
                        n = 150)
    
    contour(dens,
            add = TRUE,
            drawlabels = FALSE,
            col = "blue",
            lwd = 0.8 * facet_scale)
  }

  # -------------------------
  # Risk lines (FIXED x-range)
  # -------------------------
  x <- seq(0.001, 0.6, length.out = 200)
  
  lines(x, 0.5 / x, col = "darkred", lty = 2, lwd = 2)
  lines(x, 0.25 / x, col = "red", lty = 2, lwd = 2)
  lines(x, 0.05 / x, col = "orange", lty = 2, lwd = 2)
  
  # -------------------------
  # Legend
  # -------------------------
  legend("topright",
         legend = c("Very High Risk", "High Risk", "Moderate Risk"),
         col = c("darkred", "red", "orange"),
         lty = 2, lwd = 2 * facet_scale,
         bty = "n")
}

pdf("Plots/GRID_Mean_IDrisk.pdf",
    width = 50,
    height = 32)

facet_scale <- 1

par(
  mfrow = c(4, 5),
  mar = c(5, 4, 4, 1),
  cex = 1,
  cex.axis = 1.1 * facet_scale,
  cex.lab  = 1.3 * facet_scale,
  cex.main = 1.3 * facet_scale,
  lwd = 1 * facet_scale
)

sims <- unique(mean_across_chromosomes$Simulation)
cycles <- unique(mean_across_chromosomes$Cycle)

for (sim in sims) {
  for (cycle in cycles) {
    
    focus_run <- subset(focus_runs,
                        Simulation == sim & Cycle == cycle)
    
    bg_run <- subset(other_runs,
                     Simulation == sim & Cycle == cycle)
    
    subset_run <- subset(mean_across_chromosomes,
                          Simulation == sim & Cycle == cycle)
    
    bg_run <- bg_run[is.finite(bg_run$chrMean_FROH_2Mb) &
                       is.finite(bg_run$Mean_HET), ]
    
    focus_run <- focus_run[is.finite(focus_run$chrMean_FROH_2Mb) &
                             is.finite(focus_run$Mean_HET), ]

    subset_run <- subset_run[is.finite(subset_run$chrMean_FROH_2Mb) &
                             is.finite(subset_run$Mean_HET), ]
        
    cat(sim, cycle, "bg n =", nrow(bg_run), "\n")
    
    make_plot(bg_run, focus_run, sim, cycle)
  }
}

dev.off()

################################################################################

idrisk_results <- merged_results %>%
  mutate(Simulation = factor(Simulation, levels = c("R", "N", "W", "D")))

head(idrisk_results)

# Check for Variation
idrisk_results %>%
  group_by(Simulation, Run, Cycle) %>%
  summarise(sd_idr = sd(ID_risk_2mb, na.rm = TRUE), .groups = "drop") %>%
  summarise(min_sd = min(sd_idr), max_sd = max(sd_idr))

# Compute per-run means
idr_run_means <- idrisk_results %>%
  group_by(Simulation, Cycle, Run) %>%
  summarise(
    Run_Mean_IDR = mean(ID_risk_2mb, na.rm = TRUE),
    Run_SD_IDR   = sd(ID_risk_2mb, na.rm = TRUE),
    .groups = "drop"
  )

# Between-run variation
idr_between <- idr_run_means %>%
  group_by(Simulation, Cycle) %>%
  summarise(
    Mean_IDR = mean(Run_Mean_IDR, na.rm = TRUE),
    SD_between = sd(Run_Mean_IDR, na.rm = TRUE),
    n_runs = n(),
    SE = SD_between / sqrt(n_runs),
    
    CI_Lower = Mean_IDR - qt(0.975, df = n_runs - 1) * SE,
    CI_Upper = Mean_IDR + qt(0.975, df = n_runs - 1) * SE,
    
    .groups = "drop"
  )

# Plot
ggplot() +
  
  # each run trajectory
  geom_line(
    data = idr_run_means,
    aes(x = Cycle, y = Run_Mean_IDR, group = Run, color = Simulation),
    alpha = 0.25
  ) +
  
  # mean trajectory
  geom_line(
    data = idr_between,
    aes(x = Cycle, y = Mean_IDR, color = Simulation),
    linewidth = 1.2
  ) +
  
  # CI ribbon
  geom_ribbon(
    data = idr_between,
    aes(x = Cycle, ymin = CI_Lower, ymax = CI_Upper, fill = Simulation),
    alpha = 0.2,
    color = NA
  ) +
  
  facet_wrap(~ Simulation) +
  
  theme_minimal() +
  
  labs(
    title = "Between-Run Variation in ID Risk (2Mb)",
    x = "Cycle",
    y = "Mean ID Risk (2Mb)"
  )

ggsave("Plots/IDrisk_Run_variation.pdf", width = 8, height = 6)
ggsave("Plots/IDrisk_Run_variation.png", width = 8, height = 6)

################################################################################
## Mean IDrisk ACROSS samples and runs

mean_sim_cycle <- mean_across_chromosomes %>%
  group_by(Simulation, Cycle) %>%
  summarise(
    Mean_FROH = mean(chrMean_FROH_2Mb, na.rm = TRUE),
    Mean_HET  = mean(Mean_HET, na.rm = TRUE),
    Mean_IDR  = mean(Mean_ID_risk_2mb, na.rm = TRUE),
    .groups = "drop"
  )
head(mean_sim_cycle)

# set background point color
cycles <- sort(unique(mean_across_chromosomes$Cycle))
grey_pal <- gray(seq(0.9, 0.2, length.out = length(cycles)))
cycle_cols <- setNames(grey_pal, cycles)

# Make background cloud of points
bg_run0 <- mean_across_chromosomes %>%
  filter(Run == 0) %>%
  mutate(Grey = cycle_cols[as.character(Cycle)])

breaks <- c(0, 0.05, 0.25, 0.5, Inf)
ryPal <- colorRampPalette(c('gold1', 'darkred'))

mean_sim_cycle$Col <- ryPal(length(breaks) - 1)[
  as.numeric(cut(mean_sim_cycle$Mean_IDR,
                 breaks = breaks,
                 include.lowest = TRUE))
]

pdf("Plots/Mean_IDrisk_SIM_ONLY.pdf", width = 12, height = 9)
png("Plots/Mean_IDrisk_SIM_ONLY.png", width = 4000, height = 3000, res = 300)
png("Plots/Mean_IDrisk_SIM_ONLY.png", width = 6000, height = 5000, res = 300)

par(
  mfrow = c(2, 2),   # 4 simulations
  mar = c(5, 4, 4, 6)
)

for (sim in unique(mean_sim_cycle$Simulation)) {
  
  bg <- bg_run0[bg_run0$Simulation == sim, ]
  dat <- mean_sim_cycle[mean_sim_cycle$Simulation == sim, ]

  # --- Background (Run 0 only) ---
  plot(bg$chrMean_FROH_2Mb,
       bg$Mean_HET,
       pch = 21,
       bg  = adjustcolor(bg$Grey),
       col = NA,
       lwd = 0.3,
       cex = 1,
       xlim = c(0, 0.6),
       ylim = c(0, 2.7),
       xlab = "FROH",
       ylab = "Het/kb in non-ROH regions",
       main = paste("Simulation:", sim),
       cex.axis = 1.1,
       cex.lab = 1.3)
  
  # --- risk curves ---
  x <- seq(0.001, 0.7, length.out = 200)
  lines(x, 0.5 / x, col = "darkred", lty = 2, lwd = 2)
  lines(x, 0.25 / x, col = "red", lty = 2, lwd = 2)
  lines(x, 0.05 / x, col = "orange", lty = 2, lwd = 2)
  
  # --- main plot (ONLY mean points) ---
  dat <- dat[order(dat$Cycle), ]
  
  points(dat$Mean_FROH,
         dat$Mean_HET,
         pch = 21,
         bg  = dat$Col,
         col = "black",
         lwd = 2,
         cex = 3)

  # --- cycle labels ---
  text(dat$Mean_FROH,
       dat$Mean_HET,
       labels = dat$Cycle,
       adj = c(0.5, 0.5),
       cex = 1)
}

dev.off()
