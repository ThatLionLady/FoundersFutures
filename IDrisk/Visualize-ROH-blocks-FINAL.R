library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

m <- "25"

##### Read the PLINK output
roh_data <- read.table(paste0("ROH/m", m, ".hom"), header = TRUE)
names(roh_data) <- c('FID', 'ID', 'PHE', 'CHR', 'SNP1', 'SNP2', 'Start', 'End', 'KB', "NSNP","DENSITY", "PHOM", "PHET")
# Convert PLINK ROH length in bases from Kb to match bcftools
roh_data$Length <- roh_data$KB * 1000

chrom_names <- unique(roh_data$CHR)

################################################ PREP DATA

# Create your group column (as before):
roh_data <- roh_data %>%
  mutate(
    group = substr(ID, 1, 3),                # first 3 chars of ID
    #group = ifelse(group %in% c("SAB", "MOZ"), "MOZ", group)  # remap SAB -> MOZ
  )
head(roh_data)

# Calculate the length of each ROH in Mb
roh_data$Length_Mb <- roh_data$Length / 1e6

# Read reference genome file (chromosome sizes)
fai <- fread('References/LionRefs/LionRef.fasta.fai', header=FALSE)
names(fai) <- c('CHR', 'Chrlgth', 'Offset', 'Base', 'Width')

# Merge ROH file with chromosome lengths
roh_data <- merge(roh_data, fai[, .(CHR, Chrlgth)], by='CHR', all.x=TRUE)

# Count the number of ROH runs per individual
roh_num <- roh_data %>%
  group_by(ID) %>%
  summarise(ROH_count = n())

# Print the result
n <- length(unique(roh_data$ID))
print(n=n, roh_num)

################################################ 

# Categorize ROH by length threshold
roh_data$LengthCategory <- cut(roh_data$Length_Mb,
                                   breaks = c(0, 0.5, 1, 2, Inf),
                                   labels = c("<0.5Mb", "0.5-1Mb", "1-2Mb", ">2Mb"),
                                   right = FALSE)

####### gt1MB ####### 
# Summarize counts per individual and threshold
roh_summary <- roh_data %>%
  #Summarize counts per individual and threshold
  group_by(ID, LengthCategory) %>%
  summarise(Count = n(), .groups = 'drop')

# Display the summarized table (optional)
print(roh_summary)

# Plot the stacked bar plot
ggplot(roh_summary, aes(x = ID, y = Count, fill = factor(LengthCategory, levels = c(">2Mb", "1-2Mb", "0.5-1Mb", "<0.5Mb")))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "ROH Segment Counts by Length Category",
    x = NULL,
    y = "Number of ROH Segments",
    fill = "Length Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0("ROH/Plots/m", m, "_ROH_lengths_all.pdf"))
ggsave(paste0("ROH/Plots/m", m, "_ROH_lengths_all.png"))

####### gt500KB ####### 
# Summarize counts per individual and threshold
roh_summary <- roh_data %>%
  # Filter for ROH fragments > 500 KB (FYI: 1Mb = 1000 KB)
  filter(KB > 500) %>%
  mutate(ID = factor(ID)) %>%
  #Summarize counts per individual and threshold
  group_by(ID, LengthCategory) %>%
  summarise(Count = n(), .groups = 'drop')

# Display the summarized table (optional)
print(roh_summary)

# Plot the stacked bar plot
ggplot(roh_summary, aes(x = ID, y = Count, fill = factor(LengthCategory, levels = c(">2Mb", "1-2Mb", "0.5-1Mb", "<0.5Mb")))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "ROH Segment Counts by Length Category",
    x = NULL,
    y = "Number of ROH Segments",
    fill = "Length Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0("ROH/Plots/m", m, "_ROH_lengths_gt500KB.pdf"))
ggsave(paste0("ROH/Plots/m", m, "_ROH_lengths_gt500KB.png"))

####### gt1MB ####### 
# Summarize counts per individual and threshold
roh_summary <- roh_data %>%
  # Filter for ROH fragments > 1000 KB (1Mb)
  filter(KB > 1000) %>%
  mutate(ID = factor(ID)) %>%
  #Summarize counts per individual and threshold
  group_by(ID, LengthCategory) %>%
  summarise(Count = n(), .groups = 'drop')

# Display the summarized table (optional)
print(roh_summary)

# Plot the stacked bar plot
ggplot(roh_summary, aes(x = ID, y = Count, fill = factor(LengthCategory, levels = c(">2Mb", "1-2Mb", "0.5-1Mb", "<0.5Mb")))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "ROH Segment Counts by Length Category",
    x = NULL,
    y = "Number of ROH Segments",
    fill = "Length Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0("ROH/Plots/m", m, "_ROH_lengths_gt1MB.pdf"))
ggsave(paste0("ROH/Plots/m", m, "_ROH_lengths_gt1MB.png"))

################################################ CALCULATE FRoH

# Create results data frame
froh_results <- data.frame(Sample=character(), Group=character(), MedianChrLength=numeric(), MaxChrLength=numeric(), FROH=numeric(), stringsAsFactors=FALSE)

# Loop over each sample and calculate FROH
sample_list <- unique(roh_data$ID)
for (SAMPLE in sample_list) {
  roh_sample <- roh_data[roh_data$ID == SAMPLE, ]
  group <- unique(roh_sample$group)
  
  # Sum ROH lengths and calculate FROH for regions > 2Mb
  froh_ratio_2Mb <- sum(roh_sample$Length[roh_sample$Length > 2000000], na.rm=TRUE) / sum(fai$Chrlgth)
  
  # Sum ROH lengths and calculate FROH for regions > 2Mb
  froh_ratio_1Mb <- sum(roh_sample$Length[roh_sample$Length > 1000000], na.rm=TRUE) / sum(fai$Chrlgth)
  
  # Sum ROH lengths and calculate FROH for regions > 50Kb
  froh_ratio_50Kb <- sum(roh_sample$Length[roh_sample$Length > 50000], na.rm=TRUE) / sum(fai$Chrlgth)
  
  # Sum ROH lengths and calculate FROH regarless of size
  froh_ratio_all <- sum(roh_sample$Length, na.rm=TRUE) / sum(fai$Chrlgth)
  
  # Get median and max chromosome lengths
  median_chrlgth <- median(fai$Chrlgth)
  max_chrlgth <- max(fai$Chrlgth)
  
  # Append results to the data frame
  froh_results <- rbind(froh_results, data.frame(Sample=SAMPLE, Group=group, MedianChrLength=median_chrlgth, MaxChrLength=max_chrlgth, FROH_2Mb=froh_ratio_2Mb, FROH_1Mb=froh_ratio_1Mb, FROH_50Kb=froh_ratio_50Kb, FROH_all=froh_ratio_all))
}

summary(froh_results)
froh_results

# Write the results to FROH.txt
outfile <- "FROH.txt"

# Only write column names if the file doesn't exist yet
if (!file.exists(outfile)) {
  write.table(froh_results, file=outfile, sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
} else {
  write.table(froh_results, file=outfile, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE, append=TRUE)
}

froh_results$Sample <- factor(froh_results$Sample, levels = rev(sort(unique(froh_results$Sample))))

ggplot(froh_results, aes(y = Sample, x = FROH_1Mb, fill = Group)) +
  geom_col() +
  labs(x = expression('F'[ROH]), y = "Sample", fill = "Sample risk") +
  scale_y_discrete(limits = rev) +
  theme_minimal()

ggsave(paste0("ROH/Plots/_m", m, "_FRoH_1Mb.pdf"), width = 10, height = 5)

ggplot(froh_results, aes(y = Sample, x = FROH_2Mb, fill = Group)) +
  geom_col() +
  labs(x = expression('F'[ROH]), y = "Sample", fill = "Sample risk") +
  scale_y_discrete(limits = rev) +
  theme_minimal()

ggsave(paste0("ROH/Plots/_m", m, "_FRoH_2Mb.pdf"), width = 10, height = 5)
  
################################################ Cumulative ROH Plot

# Create bins of 10k base pairs
roh_data$LengthBin <- cut(roh_data$Length, breaks=seq(0, max(roh_data$Length), by=10000))

# Define the bin breaks manually
bin_breaks <- seq(0, max(roh_data$Length), by = 10000)

# Assign the LengthBin column based on the breaks
roh_data$LengthBin <- cut(roh_data$Length, breaks = bin_breaks, right = TRUE, include.lowest = TRUE)

# Generate the labels as '0k-10k', '10k-20k', etc.
bin_labels <- paste0((bin_breaks[-length(bin_breaks)] / 1000), "k-", (bin_breaks[-1] / 1000), "k")

# Apply the labels to the LengthBin factor
roh_data$LengthBin <- factor(roh_data$LengthBin, levels = levels(roh_data$LengthBin), labels = bin_labels)

# Calculate cumulative proportion
cumulative_data <- aggregate(Length ~ ID + LengthBin, data=roh_data, sum)
cumulative_data <- cumulative_data %>% 
  group_by(ID) %>% 
  mutate(CumProp=cumsum(Length) / sum(Length))

# Check the result
head(roh_data$LengthBin)

# Plot
ggplot(cumulative_data, aes(x=LengthBin, y=CumProp, color=ID, group=ID)) +
  geom_line() +
  labs(title="Cumulative ROH Distribution", x="ROH Length Bins", y="Cumulative Proportion") +
  theme_minimal() +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 100)]) +  # Show every 10th bin
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggsave("ROH/Plots/cumulativeROH.pdf")

### EXPLANATION:
# The plot shows the cumulative sum of ROH segments as you move from smaller ROH sizes to larger ones.
# If the line reaches 0.5 (or 50%) at a certain bin, it means that 50% of the total ROH length in the 
# dataset falls within that bin (and all the previous ones).
# Horizontal steps: Each point on the line corresponds to the total length of ROH segments in a particular 
# bin, and as you progress through bins, you accumulate the total length.

################################################ Density Plot of ROH Lengths per Sample

ggplot(roh_data, aes(x=Length, color=ID)) +
  geom_density() +
  scale_x_log10() +  # log scale if ROH lengths vary widely
  labs(title="Density of ROH Lengths across Samples", x="ROH Length (bp)", y="Density") +
  theme_minimal()

ggsave("ROH/Plots/ROH__density.pdf")

################################################################################

combined_froh_results <- read.table("FROH.txt", header = TRUE)
head(combined_froh_results)

# Plot
ggplot(combined_froh_results, aes(x = FROH_50Kb, y = FROH_1Mb, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +           # points
  #geom_line(aes(group = Group), alpha = 0.7) + # connect points per ID
  theme_minimal() +
  labs(
    title = "FROH_1Mb by Percentage per ID",
    x = "FROH_50Kb",
    y = "FROH_1Mb"
  )

ggsave(paste0("ROH/Plots/FROH_1MBv50KB.pdf"), width = 10, height = 8)

################################################################################

library(dplyr)
library(readr)

m25 <- read_table("FindingGoldilocks/imiss/m25.imiss") %>% mutate(Percentage = "25")

missing_all <- m25

head(missing_all)
summary(missing_all)

# Rename to match combined_froh_results
colnames(missing_all)[colnames(missing_all) == "INDV"] <- "Sample"
colnames(missing_all)

plot_df <- missing_all %>%
  left_join(combined_froh_results, by = c("Sample"))

sum(is.na(plot_df$FROH_1Mb))
sum(is.na(plot_df$F_MISS))

ggplot(plot_df, aes(x = F_MISS, y = FROH_1Mb, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +           # points
  #geom_line(aes(group = Group), alpha = 0.7) + # connect points per ID
  theme_minimal() +
  labs(
    title = "FROH vs Missingness per Sample",
    x = "F_MISS",
    y = "FROH_1Mb"
  )

ggsave(paste0("ROH/Plots/FROH_v_imiss.pdf"), width = 20, height = 8)
