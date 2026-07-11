################################################################################
## IDrisk

# Requires FROH.txt from Visualize-ROH-blocks.R and script calculate-heterozygosity.sh 

library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
  
################################################ PREP DATA
## Read in ROH
roh_data <- read.table("ROH/m25.hom", header = TRUE) # this version needs chromosomes changed from 1-18 to character chr (above)
names(roh_data) <- c('FID', 'ID', 'PHE', 'CHR', 'SNP1', 'SNP2', 'Start', 'End', 'KB', "NSNP","DENSITY", "PHOM", "PHET")
# Create named vector mapping numbers to names & Replace numeric chromosome values by name
chrom_names <- c("A1", "A2", "A3", "B1", "B2", "B3", "B4", "C1", "C2", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "F2", "F3")
chr_map <- setNames(chrom_names, as.character(1:18))
roh_data$CHR <- chr_map[as.character(roh_data$CHR)]
# Convert PLINK ROH length from Kb to bases to match bcftools
roh_data$Length <- roh_data$KB * 1000

## Load chromosome lengths from reference file
fai <- fread('References/LionRefs/LionRef.fasta.fai', header=FALSE)
colnames(fai) <- c('Chr', 'Chrlgth', 'Offset', 'Base', 'Width')
genome_length_kb <- sum(fai$Chrlgth) / 1000  # Convert to kilobases

################################################ Get FROH

froh_results <- read.table("FROH.txt", header = TRUE) # Calculated in Visualize-ROH-blocks.R
colnames(froh_results)
nrow(froh_results)
froh_results

################################################ Get HET-PER_KB

# if below has already been run (i.e."heterozygosity_per_kb_all_samples.csv" exists), 
# skip section below and read in data
het_results <- read.csv("heterozygosity_per_kb_all_samples.csv", header = TRUE)

################################################ PREP ANALYSIS

# Create BED directory for output
if (!dir.exists("IDrisk")) {
  dir.create("IDrisk")
}
if (!dir.exists("IDrisk/BED")) {
  dir.create("IDrisk/BED")
}
if (!dir.exists("IDrisk/HET")) {
  dir.create("IDrisk/HET")
}

################################################ CREATE BED FILE FOR EACH SAMPLE
sample_list <- unique(roh_data$ID)
total_samples <- length(sample_list)

sample_counter <- 1
for (SAMPLE in sample_list) {
  cat(sprintf("        Processing Sample %d of %d: %s\n", sample_counter, total_samples, SAMPLE))
  sample_counter <- sample_counter + 1
  
  roh_sample <- roh_data[roh_data$ID == SAMPLE, ]

  ################################################ Create BED File
  
  roh_sample <- roh_sample[roh_sample$Length > 2000000, ]
  
  # Convert to BED format
  setDT(roh_sample)
  bed_df <- roh_sample[, .(CHR, Start, End)]
  setnames(bed_df, c("CHR", "Start", "End"), c("#CHR", "START", "END"))
  
  # Save BED file
  roh_bed_file <- paste0("IDrisk/BED/roh_", SAMPLE, ".bed")
  fwrite(bed_df, file = roh_bed_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

################################################ CALCULATE nonROH HETEROZYGOSITY

## Run script to calculated nonROH-heterozygosity
system("bash Scripts/calculate-heterozygosity.sh parallel LASTfc.split.setGT-2fdp40.fixed.4pixy-m25idp800.VO.vcf.gz")

################################################ CALCULATE nonROH HETEROZYGOSITY PER KB

# Initialize an empty data frame to store all results
het_results <- data.frame()

########## VCFtools OUTPUT
# Loop through each sample and calculate heterozygosity per kb outside ROH
for (SAMPLE in sample_list) {
  # Define the input file paths
  het_file <- paste0("IDrisk/HET/noROH_indiv_", SAMPLE, ".het")
  roh_bed_file <- paste0("IDrisk/BED/roh_", SAMPLE, ".bed")
  
  # Read heterozygosity file
  het_df <- read.table(het_file, header = TRUE)
  het_df <- het_df[het_df$INDV == SAMPLE, ]
  
  # Calculate total non-ROH length
  roh_bed_df <- read.table(roh_bed_file, header = FALSE, sep = "\t")
  # Calculate non-ROH genome length in kb
  non_roh_length_kb <- genome_length_kb - sum(roh_bed_df$V3 - roh_bed_df$V2) / 1000
  
  # Compute heterozygosity per kb
  het_df$HET_PER_KB <- (het_df$N_SITES - het_df$O.HOM.) / non_roh_length_kb
  
  # Add a column for the sample ID
  het_df$Sample <- SAMPLE
  het_df$non_roh_length_kb <- non_roh_length_kb
  
  # Append the results to the all_results data frame
  het_results <- rbind(het_results, het_df)
}

# Define the output file path
output_file <- "heterozygosity_per_kb_all_samples.csv"

# Save all results to a single CSV file
write.table(het_results, file = output_file, sep = ",", quote = FALSE, row.names = FALSE)

################################################################################
# SKIP TO HERE...

# Create the stacked bar plot
ggplot(het_results, aes(x = Sample, y = HET_PER_KB)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Het/kb in non−ROH regions",
       x = "Sample ID",
       y = "Het/kb") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

################################################ PLOT IDrisk

# Comparative species data from Chris
idrisk_data <- read.csv("Scripts/ROH_het_by_species.csv")
names(idrisk_data) <- c("Sample", "HET_PER_KB", "FROH_2Mb", "ID_risk_2mb")
idrisk_data$Sample

# Merge FROH and heterozygosity results from FROH.R and nonROH-het-per-kb.R
merged_results <- merge(froh_results[, c("Sample", "FROH_2Mb")], het_results[, c("Sample", "HET_PER_KB")], by = "Sample")

# Calculate ID risk
merged_results <- merged_results %>%
  mutate(ID_risk_2mb = HET_PER_KB * FROH_2Mb)

# Combine with idrisk_data, keeping only relevant columns for plotting
data <- rbind(idrisk_data, merged_results)

# Ensure ID_risk_2mb is numeric in the merged data
data$ID_risk_2mb <- as.numeric(data$ID_risk_2mb)

# Define the list of samples to exclude
exclude_points <- c(
  #"Florida panther",
  "Chatham island black robin",
  #"Isle Royale wolf",   
  #"Los Angeles mountain lion",
  "Pacific pocket mouse",
  "Kakapo (stewart island)",
  "California condor",
  #"Cape Peninsula caracals",
  #"Minnesota wolf",
  "Mountain gorilla",
  "Seychelles magpie robin",
  "Orange bellied parrot",
  #"Channel island fox",
  "Northern white rhino",
  "Gulf of California Fin Whales",
  "Southern resident killer whale",
  "Wrangel Island Woolly Mammoth",
  "'Akikiki",
  "Alpine ibex",
  "Isle Royale moose",
  #"Ethiopian wolf",
  "Minnesota moose",
  "Vaquita",
  "Poouli",
  "Whooping crane"
)

# Exclude specific samples from idrisk_data
data <- data %>% filter(!Sample %in% exclude_points)

################################################ MY PLOT BASED ON Chris's

pdf("Plots/plot_ID_risk_2Mb.pdf", width = 8, height = 6)

ryPal <- colorRampPalette(c('gold1','dark red'))
data$Col <- ryPal(20)[as.numeric(cut(data$ID_risk_2mb,breaks = 20))]
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))

plot(data$FROH_2Mb, data$HET_PER_KB, col = data$Col, pch=19, xlim = c(0, 0.6), ylim = c(0, 2.7), 
     ylab = "Het/kb in non-ROH regions", xlab =  "FROH", cex.axis=1.1, cex.lab=1.3)
text(data$FROH_2Mb, data$HET_PER_KB, labels = data$Sample, cex=0.5, adj=c(0, 0.5), srt=-20) # Adjust cex for text size and adj for position
x <- seq(0, 0.7, length.out = 100)
vhigh_risk <- 0.5 / x
high_risk <- 0.25 / x
moderate_risk <- 0.05/x
lines(x,vhigh_risk, col = data$Col[2])
lines(x,high_risk, col = data$Col[7])
lines(x,moderate_risk, col = data$Col[15])

colors <- ryPal(20)
par(mar = c(5, 4, 4, 6)) # Adjust margins to make space for the legend
image(1, seq(0, 0.9, length.out = 100), t(seq(0, 0.9, length.out = 100)), col = colors, axes = FALSE, xlab = "", ylab = "")
axis(4, at = seq(0, 0.9, by = .1), labels = seq(0, 0.9, by = 0.1), las = 1)
mtext("IDrisk = FROH * non-ROH Heterozygosity", side = 4, line = 2.5)

dev.off()

################################################ MY GGPLOT PLOT

# Create a data frame for risk lines
risk_lines <- data.frame(
  x = seq(0.01, 0.7, length.out = 100),
  vhigh_risk = 0.5 / seq(0.01, 0.7, length.out = 100),
  high_risk = 0.25 / seq(0.01, 0.7, length.out = 100),
  moderate_risk = 0.05 / seq(0.01, 0.7, length.out = 100)
)

# Plotting with ggplot2
ggplot(data, aes(x = FROH_2Mb, y = HET_PER_KB, color = ID_risk_2mb)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "gold1", high = "darkred") +
  labs(x = expression('F'[ROH]), y = "Het/kb in non-ROH regions", color = "ID risk") +
  theme_minimal() +
  geom_text(aes(label = paste0(Sample, "\n", round(ID_risk_2mb, 2))), hjust = 0, vjust = -1, size = 3, check_overlap = TRUE) +
  xlim(0, 0.6) + 
  ylim(0, 2.7) +
  geom_line(data = risk_lines, aes(x = x, y = vhigh_risk), color = "red") +
  geom_line(data = risk_lines, aes(x = x, y = high_risk), color = "orange") +
  geom_line(data = risk_lines, aes(x = x, y = moderate_risk), color = "yellow")

# Highlight specific sample

ggplot(data, aes(x = FROH_2Mb, y = HET_PER_KB, color = ID_risk_2mb)) +
  scale_color_gradient(low = "gold1", high = "darkred") +
  labs(x = expression('F'[ROH]), y = "Het/kb in non-ROH regions", color = "ID risk") +
  theme_minimal() +
  geom_text(aes(label = paste0(Sample, "\n", round(ID_risk_2mb, 2))), hjust = 0, vjust = -1, size = 3, check_overlap = TRUE) +
  xlim(0, 0.6) + 
  ylim(0, 2.7) +
  geom_line(data = risk_lines, aes(x = x, y = vhigh_risk), color = "red") +
  geom_line(data = risk_lines, aes(x = x, y = high_risk), color = "orange") +
  geom_line(data = risk_lines, aes(x = x, y = moderate_risk), color = "yellow") +
  geom_point(size = 3) +
  geom_point(data = subset(data, Sample == "TSW-M002"), aes(x = FROH_2Mb, y = HET_PER_KB), 
             size = 3, color = "black", stroke = 1)

#################################################################################
## Single Value for ZDC Lions

library(dplyr)
library(stringr)
library(ggrepel)

zdc_summary <- data %>%
  mutate(
    is_ZDC = str_detect(Sample, "^[A-Z]{3}-"),
    Pop = ifelse(is_ZDC, str_extract(Sample, "^[A-Z]{3}"), Sample)
  ) %>%
  filter(is_ZDC) %>%
  group_by(Pop) %>%
  reframe(
    Sample = Pop[1],
    HET_PER_KB  = mean(HET_PER_KB, na.rm = TRUE),
    FROH_2Mb    = mean(FROH_2Mb, na.rm = TRUE),
    ID_risk_2mb = mean(ID_risk_2mb, na.rm = TRUE),
    Col = first(Col)
  )

final_data <- bind_rows(idrisk_data, zdc_summary)

# Define the list of samples to exclude
exclude_points <- c(
  #"Florida panther",
  "Chatham island black robin",
  #"Isle Royale wolf",   
  #"Los Angeles mountain lion",
  "Pacific pocket mouse",
  "Kakapo (stewart island)",
  "California condor",
  #"Cape Peninsula caracals",
  #"Minnesota wolf",
  "Mountain gorilla",
  "Seychelles magpie robin",
  "Orange bellied parrot",
  #"Channel island fox",
  "Northern white rhino",
  "Gulf of California Fin Whales",
  "Southern resident killer whale",
  "Wrangel Island Woolly Mammoth",
  "'Akikiki",
  "Alpine ibex",
  "Isle Royale moose",
  #"Ethiopian wolf",
  "Minnesota moose",
  "Vaquita",
  "Poouli",
  "Whooping crane"
)

# Exclude specific samples from idrisk_data
final_data <- final_data %>% filter(!Sample %in% exclude_points)

colnames(final_data)
final_data

# Create a data frame for risk lines
risk_lines <- data.frame(
  x = seq(0.01, 0.7, length.out = 100),
  vhigh_risk = 0.5 / seq(0.01, 0.7, length.out = 100),
  high_risk = 0.25 / seq(0.01, 0.7, length.out = 100),
  moderate_risk = 0.05 / seq(0.01, 0.7, length.out = 100)
)

# Plotting with ggplot2
ggplot(final_data, aes(x = FROH_2Mb, y = HET_PER_KB, fill = ID_risk_2mb)) +
  geom_line(data = risk_lines, aes(x = x, y = vhigh_risk), color = "red", inherit.aes = FALSE) +
  geom_line(data = risk_lines, aes(x = x, y = high_risk), color = "orange", inherit.aes = FALSE) +
  geom_line(data = risk_lines, aes(x = x, y = moderate_risk), color = "yellow", inherit.aes = FALSE) +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.8) +
  scale_fill_gradient(low = "gold1", high = "darkred") +
  labs(x = expression('F'[ROH]), y = "Het/kb in non-ROH regions", fill = "ID risk") +
  theme_minimal() +
  #geom_text(aes(label = paste0(Sample, " (", round(ID_risk_2mb, 3), ")")), color = "black", hjust = 0.5, vjust = -1, size = 4) +
  geom_text_repel(
    aes(label = paste0(Sample, " (", round(ID_risk_2mb, 3), ")")),
    color = "black",
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  xlim(0, 0.6) + 
  ylim(0, 2) +
  theme_grey() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")) 

ggsave("Plots/IDrisk_bySource.pdf", width = 8, height = 8)
ggsave("Plots/IDrisk_bySource.png", width = 8, height = 8)
