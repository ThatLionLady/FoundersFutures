library(tidyverse)

setwd("ZDC/SLiM/")

# Collect files
eig_files <- list.files("PCAs", pattern = "^A1_.*\\.eigenvec$", full.names = TRUE)
length(eig_files)

# Read in files
library(tidyverse)

files <- list.files("PCAs", pattern = "^A1_.*\\.eigenvec$", full.names = TRUE)

eig_complete <- map_dfr(eig_files, function(f) {
  
  base <- basename(f)
  
  # extract components directly
  m <- str_match(base, "^A1_([WNDR])([0-9]+)_([0-9]+)\\.eigenvec$")
  
  tibble(
    SimType = m[2],   # D or R
    Run     = as.integer(m[3]),
    Cycle   = as.integer(m[4])
  ) %>%
    bind_cols(
      read_tsv(f, col_names = TRUE)
    )
})
head(eig_complete)

# -----------------------------
# Add input PCA files
# -----------------------------

# Reintroduced -> R Cycle 1
reintro <- read_tsv(
  "ZDC/SLiM/Input/PCAs/Reintroduced.vcf.eigenvec",
  col_names = TRUE
) %>%
  mutate(
    SimType = factor("R", levels = levels(eig_complete$SimType)),
    Run     = factor(0, levels = levels(eig_complete$Run)),
    Cycle   = 1
  )

# TrueFounders -> N/W/D Cycle 1
founders_raw <- read_tsv(
  "ZDC/SLiM/Input/PCAs/TrueFounders.eigenvec",
  col_names = TRUE
)

founders <- map_dfr(
  c("N", "W", "D"),
  ~ founders_raw %>%
    mutate(
      SimType = factor(.x, levels = levels(eig_complete$SimType)),
      Run     = factor(0, levels = levels(eig_complete$Run)),
      Cycle   = 1
    )
)

# Combine everything
eig_complete <- bind_rows(
  eig_complete,
  reintro,
  founders
)

# what did we successfully parse?
eig_complete %>%
  count(SimType, Run, Cycle) %>%
  arrange(SimType, Run, Cycle) %>%
  print(n = 160)

eig_complete <- eig_complete %>%
  mutate(
    SimType = factor(SimType, levels = c("R", "N", "W", "D")),
    Run = factor(Run, levels = 0:9)
  )

# Reduce background points
eig_plot <- eig_complete %>%
  group_by(SimType, Cycle) %>%
  group_modify(~ {
    
    bg <- filter(.x, Run != 0)
    fg <- filter(.x, Run == 0)
    
    # keep 25% of background points
    bg_s <- slice_sample(bg, prop = 0.25)
    
    bind_rows(bg_s, fg)
    
  }) %>%
  ungroup()

######################
# Overlapping run plot (single run focus w background dots)
eig_plot$plot_group <- ifelse(
  eig_plot$Run == 9,
  as.character(eig_plot$Cluster),
  "grey85"
)

# PC1 - PC2
ggplot(eig_plot, aes(x = PC1, y = PC2)) +
  
  # Background points
  geom_point(
    data = subset(eig_plot, Run != 0),
    color = "grey80",
    alpha = 0.3
  ) +
  
  # Highlighted Run0 points
  geom_point(
    data = subset(eig_plot, Run == 0),
    aes(color = Cluster),
    alpha = 0.8
  ) +
  
  facet_grid(SimType ~ Cycle) +
  theme_minimal()

ggsave("Plots/allPCAs_Overlay_PC1PC2.pdf", width = 20, height = 5)
ggsave("Plots/allPCAs_Overlay_PC1PC2.png", width = 20, height = 5)

# PC1 - PC3
ggplot(eig_plot, aes(x = PC1, y = PC3)) +
  
  # Background points
  geom_point(
    data = subset(eig_plot, Run != 0),
    color = "grey80",
    alpha = 0.3
  ) +
  
  # Highlighted Run0 points
  geom_point(
    data = subset(eig_plot, Run == 0),
    aes(color = Cluster),
    alpha = 0.8
  ) +
  
  facet_grid(SimType ~ Cycle) +
  theme_minimal()

ggsave("Plots/allPCAs_Overlay_PC1PC3.pdf", width = 20, height = 5)
ggsave("Plots/allPCAs_Overlay_PC1PC3.png", width = 20, height = 5)


######################
# Overlapping run plot (no focus)
ggplot(
  eig_plot,
  aes(x = PC1, y = PC3, color = factor(Run))) +
  geom_point(alpha = 0.3) +
  facet_grid(SimType ~ Cycle) +
  theme_minimal()

ggsave("Plots/allPCAs_Overlay.pdf", width = 20, height = 5)
ggsave("Plots/allPCAs_Overlay.png", width = 20, height = 5)

######################
# 3D plot with run focus

library(tidyverse)

# Example subset (TEST 3d plotting with single sim/cycle)
plot_df <- eig_plot %>%
  filter(SimType == "W",
         Cycle == 1)

# ----- Perspective projection -----

project_3d <- function(x, y, z,
                       angle_x = pi/6,
                       angle_y = pi/4,
                       z_scale = 1) {
  
  x2 <- x*cos(angle_y) + y*sin(angle_y)
  
  y2 <- z*z_scale +
    x*sin(angle_x)*sin(angle_y) -
    y*sin(angle_x)*cos(angle_y)
  
  tibble(x2, y2)
}

proj <- project_3d(
  plot_df$PC1,
  plot_df$PC2,
  plot_df$PC3
)

plot_df <- bind_cols(plot_df, proj)

# ----- Axes -----

axis_len <- max(abs(c(
  plot_df$PC1,
  plot_df$PC2,
  plot_df$PC3
)))

# positive + negative ends
axes <- tribble(
  ~x, ~y, ~z, ~xend, ~yend, ~zend, ~label, ~labx, ~laby, ~labz,
  
  -axis_len, 0, 0,  axis_len, 0, 0, "PC1",  axis_len, 0, 0,
  0,-axis_len, 0,  0, axis_len, 0, "PC2",  0, axis_len, 0,
  0, 0,-axis_len,  0, 0, axis_len, "PC3",  0, 0, axis_len
)

# project line starts
start_proj <- project_3d(
  axes$x,
  axes$y,
  axes$z
)

# project line ends
end_proj <- project_3d(
  axes$xend,
  axes$yend,
  axes$zend
)

# project label positions
label_proj <- project_3d(
  axes$labx,
  axes$laby,
  axes$labz
)

axes <- bind_cols(
  axes,
  start_proj |> rename(x_start = x2, y_start = y2),
  end_proj   |> rename(x_end   = x2, y_end   = y2),
  label_proj |> rename(x_lab   = x2, y_lab   = y2)
)

# ----- Plot -----

ggplot(plot_df, aes(x2, y2)) +
  
  # axes
  geom_segment(
    data = axes,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end
    ),
    linewidth = 0.5,
    color = "grey50"
  ) +
  
  # labels
  geom_text(
    data = axes,
    aes(
      x = x_lab,
      y = y_lab,
      label = label
    ),
    size = 4,
    color = "grey30"
  ) +
  
  # points
  geom_point(
    data = subset(plot_df, Run != 0),
    color = "grey80",
    alpha = 0.3,
    size = 2
  ) +
  
  geom_point(
    data = subset(plot_df, Run == 0),
    aes(color = Cluster),
    size = 2.5
  ) +
  
  coord_equal() +
  theme_void()

###################
## Full Facet Grid (all Sim/Cycles)

library(tidyverse)

# ----- Projection function -----

project_3d <- function(x, y, z,
                       angle_x = pi/7,
                       angle_y = pi/4,
                       z_scale = 1) {
  
  x2 <- x*cos(angle_y) + y*sin(angle_y)
  
  y2 <- z*z_scale +
    x*sin(angle_x)*sin(angle_y) -
    y*sin(angle_x)*cos(angle_y)
  
  tibble(x2, y2)
}

# ----- Project points -----

proj <- project_3d(
  eig_plot$PC1,
  eig_plot$PC2,
  eig_plot$PC3
)

plot_df <- bind_cols(eig_plot, proj)

# ----- Build axes for every facet -----

axis_len <- max(abs(c(
  eig_plot$PC1,
  eig_plot$PC2,
  eig_plot$PC3
)))

facet_combos <- expand_grid(
  SimType = unique(plot_df$SimType),
  Cycle = unique(plot_df$Cycle)
)

axes_base <- tribble(
  ~x, ~y, ~z, ~xend, ~yend, ~zend, ~label, ~labx, ~laby, ~labz,
  
  -axis_len, 0, 0,  axis_len, 0, 0, "PC1",  axis_len, 0, 0,
  0,-axis_len, 0,  0, axis_len, 0, "PC2",  0, axis_len, 0,
  0, 0,-axis_len,  0, 0, axis_len, "PC3",  0, 0, axis_len
)

axes <- facet_combos %>%
  crossing(axes_base)

# project axis starts
start_proj <- project_3d(
  axes$x,
  axes$y,
  axes$z
)

# project axis ends
end_proj <- project_3d(
  axes$xend,
  axes$yend,
  axes$zend
)

# project labels
label_proj <- project_3d(
  axes$labx,
  axes$laby,
  axes$labz
)

axes <- bind_cols(
  axes,
  start_proj |> rename(x_start = x2, y_start = y2),
  end_proj   |> rename(x_end   = x2, y_end   = y2),
  label_proj |> rename(x_lab   = x2, y_lab   = y2)
)

# ----- Plot -----

ggplot(plot_df, aes(x2, y2)) +
  
  # axes
  geom_segment(
    data = axes,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end
    ),
    linewidth = 0.5,
    color = "grey50",
    inherit.aes = FALSE
  ) +
  
  # axis labels
  geom_text(
    data = axes,
    aes(
      x = x_lab,
      y = y_lab,
      label = label
    ),
    size = 2,
    color = "grey30",
    inherit.aes = FALSE
  ) +
  
  # background points
  geom_point(
    data = subset(plot_df, Run != 0),
    color = "grey85",
    shape = 16,
    alpha = 0.4,
    size = 0.5
  ) +
  
  # highlighted points
  geom_point(
    data = subset(plot_df, Run == 0),
    aes(color = Cluster),
    shape = 16,
    size = 1
  ) +
  
  facet_grid(SimType ~ Cycle) +
  
  coord_equal() +
  
  theme_void() +
  
  theme(
    #legend.position = "none",
    strip.text = element_text(size = 12),
    panel.spacing = unit(1, "lines")
  )

ggsave("Plots/3D_PCAs.png", width = 20, height = 5)
ggsave("Plots/3D_PCAs-legend.pdf", width = 20, height = 5, device = cairo_pdf)

######################
# Single Run Focus

for (run in 0:9) { 

  ggplot(
    eig_plot %>% 
      filter(Run == run),
         aes(x = PC1, y = PC2, color = factor(Cluster))) +
    geom_point(alpha = 0.3) +
    facet_grid(SimType ~ Cycle) +
    theme_minimal()
  
  ggsave(paste0("Plots/PCA_Run", run, "_Overlay.pdf"), width = 20, height = 5)
  ggsave(paste0("Plots/PCA_Run", run, "_Overlay.png"), width = 20, height = 5)

}

### --- CLUSTER COUNT --- ###

library(dplyr)

cluster_summary <- eig_plot %>%
  group_by(SimType, Cycle, Run) %>%
  summarise(n_clusters = n_distinct(Cluster), .groups = "drop")

cluster_summary <- cluster_summary %>%
  mutate(
    SimType = factor(SimType, levels = c("R", "N", "W", "D")),
    Run = factor(Run, levels = 0:9)
  )

cluster_summary <- cluster_summary %>%
  mutate(
    x_order = interaction(SimType, Run, sep = "")
  ) %>%
  arrange(SimType, Run) %>%
  mutate(x_order = factor(x_order, levels = unique(x_order)))

write_csv(cluster_summary, "cluster_summary.csv")

library(ggplot2)

# Box plot
ggplot(cluster_summary,
       aes(x = factor(Run), y = n_clusters, fill = SimType)) +
  geom_col() +
  geom_smooth(aes(group = SimType),
              color = "black",
              method = "loess",
              se = FALSE,
              linewidth = 0.5) +
  facet_grid(SimType~ Cycle) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust = 1)
  ) +
  xlab("Run") +
  ylab("Number of Clusters")

ggsave("Plots/allPCAs_box_plot.pdf", width = 20, height = 5)
ggsave("Plots/allPCAs_box_plot.png", width = 20, height = 5)

# Line plot
ggplot(cluster_summary,
       aes(x = factor(Run),
           y = n_clusters,
           color = SimType,
           group = SimType)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_grid(SimType ~ Cycle) +
  theme_bw() +
  xlab("Run") +
  ylab("Number of Clusters")

ggsave("Plots/allPCAs_line_plot.pdf", width = 20, height = 5)
ggsave("Plots/allPCAs_line_plot.png", width = 20, height = 5)