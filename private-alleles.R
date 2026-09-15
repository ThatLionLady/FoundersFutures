library(poppr)
library(ggplot2)
library(vcfR)

############################################################################################
## Get Data

#remotes::install_github("dinmatias/reconproGS")
library(reconproGS)
library(adegenet)

m <- "m00"

# Load VARIANT ONLY VCF as genind
vcf <- read.vcfR(paste0("LASTfc.split.setGT-2fdp40.fixed.4pixy-", m, "idp800.VO.vcf.gz"), verbose = FALSE)
vcfGenind_ini <- vcfR2genind(vcf)

# Count loci in genind
nLoc(vcfGenind_ini)

vcfGenind <- vcfGenind_ini

## Add Pop info to genind NEW WAY : USE THIS!!!! ##

# Extract population code from individual names (first 3 letters)
pop_codes <- substr(indNames(vcfGenind), 1, 3)

# Assign to genind
vcfGenind@pop <- as.factor(pop_codes)

# Verify alignment
verification_df <- data.frame(
  Individual = indNames(vcfGenind),
  Population = pop(vcfGenind)
)

# View the data frame to check that ID and pop match
print(verification_df)

############################################################################################
# Print number of private alleles per POPULATION
# Creates transposed_matrix

# Obtain a list of private alleles for each population
private_alleles_matrix <- private_alleles(vcfGenind)
# Print the first 5 rows and 5 columns of the matrix
head(private_alleles_matrix, 6)[, 1:5]

# Transpose Matrix
transposed_matrix <- t(private_alleles_matrix)
head(transposed_matrix)

dim(transposed_matrix)
summary(transposed_matrix)
colSums(transposed_matrix, na.rm = TRUE)

###### CHECK DISTRIBUTION

# Sum allele counts across all populations for each SNP
total_allele_counts <- rowSums(transposed_matrix)

# Convert to a dataframe for plotting
total_allele_counts_df <- data.frame(
  SNP = names(total_allele_counts),
  Allele_Count = total_allele_counts
)

# Plot distribution of total allele counts
ggplot(total_allele_counts_df, aes(x = Allele_Count)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of Allele Counts Across All Individuals",
    x = "Total Allele Count",
    y = "Number of SNPs"
  ) +
  theme_minimal()

ggsave(paste0("Plots/PA/", m, "_distribution_of_allele_counts.pdf"), width = 10, height = 8)
ggsave(paste0("Plots/PA/", m, "_distribution_of_allele_counts.png"), width = 10, height = 8)

# Why do these zeros exist? 
zero_rows <- rowSums(transposed_matrix) == 0
rownames(transposed_matrix)[which(zero_rows)]
transposed_matrix[zero_rows, , drop = FALSE]

sum(zero_rows)  # number of alleles with 0 counts
# [1] 0
sum(transposed_matrix)

# Actual SNP count (bcftools view -H | wc -l)

############################################################################################
## Fixed Alleles (YAY IT'S PERFECT!)
# Creates filtered_matrix REQUIRED

# Load the required package
library(adegenet)

#### Get Population Counts ####

# Get the population information for each individual
populations <- pop(vcfGenind)

# Count the number of individuals in each population & 
# convert to a data frame for easier viewing
population_counts_df <- as.data.frame(table(populations))
colnames(population_counts_df) <- c("Population", "Individual_Count")

# View the result
print(population_counts_df)

# Convert population_counts_df to a named vector for easier comparison
population_counts <- setNames(population_counts_df$Individual_Count, population_counts_df$Population)
population_diploid_counts <- population_counts * 2  # diploid correction
population_diploid_counts

#### Collect Fixed Alleles ####

# Filter rows in transposed_matrix where values match the population count * 2 (i.e. FIXED alleles)

filtered_matrix <- as.data.frame(transposed_matrix)

filtered_matrix <- filtered_matrix[
  (filtered_matrix$KAR == population_diploid_counts["KAR"] |
     filtered_matrix$KHA == population_diploid_counts["KHA"] |
     filtered_matrix$MAK == population_diploid_counts["MAK"] |
     filtered_matrix$MOZ == population_diploid_counts["MOZ"] |
     filtered_matrix$SAB == population_diploid_counts["SAB"] |
     filtered_matrix$TEM == population_diploid_counts["TEM"] |
     filtered_matrix$TSW == population_diploid_counts["TSW"]),
]

# View the filtered matrix
head(filtered_matrix)
summary(filtered_matrix)
sum(colSums(filtered_matrix != 0)) # sanity check
nrow(filtered_matrix) # should equal above

## Extra sanity check...
# Assume filtered_matrix is your matrix of SNPs x populations
# Convert to data frame
filtered_df <- as.data.frame(filtered_matrix)

# Create an empty list to store results per population
allele_count_freq <- list()

library(dplyr)

for(pop in colnames(filtered_df)) {
  # Count how many SNPs have each allele count
  freq_table <- filtered_df %>%
    count(Allele_Count = !!sym(pop)) %>%   # count occurrences of each allele count
    arrange(Allele_Count)
  
  # Add a column for the population name
  freq_table$Population <- pop
  
  allele_count_freq[[pop]] <- freq_table
}

# Combine all populations into one data frame
allele_count_freq_df <- bind_rows(allele_count_freq)

# Optional: remove zeros if you don't want them
allele_count_freq_df_no0 <- allele_count_freq_df %>% filter(Allele_Count > 0)

# View
print(allele_count_freq_df_no0) # <- all of these should equal the population_diploid_counts

###### THIS WILL BE THE END OF ME! ######
### Calling it good... Let's do this...

## Make a data frame of allele type/condition counts 
# Count fixed private alleles per population
fixed_private_alleles_count <- colSums(filtered_matrix != 0)
print(fixed_private_alleles_count)

# Count total private alleles per population
total_private_alleles_count <- colSums(transposed_matrix != 0)
print(total_private_alleles_count)

# Combine results into a data frame
allele_summary <- data.frame(
  Population = colnames(transposed_matrix),
  Fixed_Private_Count = fixed_private_alleles_count,
  Total_Private_Count = total_private_alleles_count
)

# Add non-fixed private alleles (Total - Fixed)
allele_summary$Non_Fixed_Count <- allele_summary$Total_Private_Count - allele_summary$Fixed_Private_Count

# If there are populations with ZERO fixed alleles
allele_summary <- allele_summary %>%
  mutate(
    Fixed_Private_Count = as.numeric(ifelse(is.na(Fixed_Private_Count), 0, Fixed_Private_Count)),
    Total_Private_Count = as.numeric(ifelse(is.na(Total_Private_Count), 0, Total_Private_Count))
  )

# View the summary
print(allele_summary)

write.table(allele_summary, paste0(m, "_private-alleles-by-source.txt"))

#### PLOT

# Reshape the data to long format for plotting
population_counts_long <- allele_summary %>%
  pivot_longer(
    cols = c(Fixed_Private_Count,Non_Fixed_Count),
    names_to = "Condition",
    values_to = "Count"
  )

print(population_counts_long)

# Plot the data
ggplot(population_counts_long, aes(fill = Condition, y = Count, x = Population)) + 
  geom_bar(position = "stack", stat = "identity") +
  labs(
    title = "Private Alleles by Population",
    x = "Population",
    y = "Number of Alleles"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("Plots/PA/", m, "_Private_Alleles_by_Source.pdf"), width = 10, height = 8)
ggsave(paste0("Plots/PA/", m, "_Private_Alleles_by_Source.png"), width = 10, height = 8)

# Filter to include only "Fixed" counts
fixed_counts <- population_counts_long %>%
  filter(Condition == "Fixed_Private_Count")

# Plot the data
ggplot(fixed_counts, aes(fill = Condition, y = Count, x = Population)) + 
  geom_bar(position = "stack", stat = "identity") +
  labs(
    title = "Fixed Alleles by Population",
    x = "Population",
    y = "Number of Fixed Alleles"
  ) +
  theme_minimal() +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("Plots/PA/", m, "_Fixed_Alleles_by_Source.pdf"), width = 10, height = 8)
ggsave(paste0("Plots/PA/", m, "_Fixed_Alleles_by_Source.png"), width = 10, height = 8)


# ===============================
# FIXED PRIVATE ALLELES PER ORIGIN
# ===============================

# Origin membership table
origins_df <- read.table("origins-pops.txt", header = TRUE)
origins_df$Kgalagadi <- origins_df$Kgalagahi
summary(origins_df)

### KRUGER

# Pick the origin of interest
origin_name <- "Kruger"
vcfGenind_KR <- vcfGenind

# Extract the genind individual names
genind_ids <- indNames(vcfGenind_KR)

# Reorder origins_df to match genind IDs
origins_ordered <- origins_df[match(genind_ids, origins_df$ID), ]

# Safety check: ensure ordering worked
stopifnot(all(origins_ordered$ID == genind_ids))

# Create a new population factor:
# Individuals that are TRUE for this origin get the origin name
# All others are "Other"
KR_pop_codes <- ifelse(origins_ordered[[origin_name]], origin_name, "Other")

# Assign it to the genind object
vcfGenind_KR@pop <- as.factor(KR_pop_codes)

# Check
table(pop(vcfGenind_KR))

# Verify alignment
verification_df <- data.frame(
  Individual = indNames(vcfGenind_KR),
  Population = pop(vcfGenind_KR)
)

# View the data frame to check that ID and pop match
print(verification_df)

# Obtain a list of private alleles for each population
KR_private_alleles_matrix <- private_alleles(vcfGenind_KR)
# Print the first 5 rows and 5 columns of the matrix
head(KR_private_alleles_matrix, 6)[, 1:5]

# Transpose Matrix
KR_transposed_matrix <- t(KR_private_alleles_matrix)
head(KR_transposed_matrix)

dim(KR_transposed_matrix)

### ETOSHA

# Pick the origin of interest
origin_name <- "Etosha"
vcfGenind_ET <- vcfGenind

# Extract the genind individual names
genind_ids <- indNames(vcfGenind_ET)

# Reorder origins_df to match genind IDs
origins_ordered <- origins_df[match(genind_ids, origins_df$ID), ]

# Safety check: ensure ordering worked
stopifnot(all(origins_ordered$ID == genind_ids))

# Create a new population factor:
# Individuals that are TRUE for this origin get the origin name
# All others are "Other"
ET_pop_codes <- ifelse(origins_ordered[[origin_name]], origin_name, "Other")

# Assign it to the genind object
vcfGenind_ET@pop <- as.factor(ET_pop_codes)

# Check
table(pop(vcfGenind_ET))

# Verify alignment
verification_df <- data.frame(
  Individual = indNames(vcfGenind_ET),
  Population = pop(vcfGenind_ET)
)

# View the data frame to check that ID and pop match
print(verification_df)

# Obtain a list of private alleles for each population
ET_private_alleles_matrix <- private_alleles(vcfGenind_ET)
# Print the first 5 rows and 5 columns of the matrix
head(ET_private_alleles_matrix, 6)[, 1:5]

# Transpose Matrix
ET_transposed_matrix <- t(ET_private_alleles_matrix)
head(ET_transposed_matrix)

dim(ET_transposed_matrix)

### Kgalagadi

# Pick the origin of interest
origin_name <- "Kgalagadi" # spelled wrong (whoops!) CHANGE BEFORE PUBLISHING!
vcfGenind_KG <- vcfGenind

# Extract the genind individual names
genind_ids <- indNames(vcfGenind_KG)

# Reorder origins_df to match genind IDs
origins_ordered <- origins_df[match(genind_ids, origins_df$ID), ]

# SafKGy check: ensure ordering worked
stopifnot(all(origins_ordered$ID == genind_ids))

# Create a new population factor:
# Individuals that are TRUE for this origin get the origin name
# All others are "Other"
KG_pop_codes <- ifelse(origins_ordered[[origin_name]], origin_name, "Other")

# Assign it to the genind object
vcfGenind_KG@pop <- as.factor(KG_pop_codes)

# Check
table(pop(vcfGenind_KG))

# Verify alignment
verification_df <- data.frame(
  Individual = indNames(vcfGenind_KG),
  Population = pop(vcfGenind_KG)
)

# View the data frame to check that ID and pop match
print(verification_df)

# Obtain a list of private alleles for each population
KG_private_alleles_matrix <- private_alleles(vcfGenind_KG)
# Print the first 5 rows and 5 columns of the matrix
head(KG_private_alleles_matrix, 6)[, 1:5]

# Transpose Matrix
KG_transposed_matrix <- t(KG_private_alleles_matrix)
head(KG_transposed_matrix)

dim(KG_transposed_matrix)

# Filter to just private alleles per origin

KR_private_only <- KR_transposed_matrix[KR_transposed_matrix[, "Kruger"] != 0, ]
ET_private_only <- ET_transposed_matrix[ET_transposed_matrix[, "Etosha"] != 0, ]
KG_private_only <- KG_transposed_matrix[KG_transposed_matrix[, "Kgalagadi"] != 0, ]

# Combine lists and remove duplicates

# Combine the private-only matrices into data frames with SNP identifiers
KR_df <- data.frame(SNP = rownames(KR_private_only), Kruger = KR_private_only[, "Kruger"])
ET_df <- data.frame(SNP = rownames(ET_private_only), Etosha = ET_private_only[, "Etosha"])
KG_df <- data.frame(SNP = rownames(KG_private_only), Kgalagadi = KG_private_only[, "Kgalagadi"])

# Combine them with full joins, merging on SNP
combined_origins <- KR_df %>%
  full_join(ET_df, by = "SNP") %>%
  full_join(KG_df, by = "SNP") %>%
  # Replace NAs with 0s
  mutate(across(-SNP, ~replace_na(., 0)))

# Set SNP as rownames if you prefer a "matrix-style" object
combined_origins_matrix <- as.data.frame(combined_origins)
rownames(combined_origins_matrix) <- combined_origins_matrix$SNP
combined_origins_matrix <- combined_origins_matrix[, -1]

# Check result
dim(combined_origins_matrix)
head(combined_origins_matrix)

# Delete none-private rows

# Keep only SNPs private to exactly one origin
origins_private_only_matrix <- combined_origins_matrix[rowSums(combined_origins_matrix > 0) == 1, ]

# Check dimensions and preview
dim(origins_private_only_matrix)
head(origins_private_only_matrix)
summary(origins_private_only_matrix)
sapply(origins_private_only_matrix, range)

# Are there any fixed alleles?

# Convert to tidy long format
origins_private_long <- origins_private_only_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SNP") %>%
  pivot_longer(cols = c("Kruger", "Etosha", "Kgalagadi"),
               names_to = "Origin",
               values_to = "Allele_Count") %>%
  filter(Allele_Count > 0)

# Plot histograms of allele counts per origin
ggplot(origins_private_long, aes(x = Allele_Count)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  facet_wrap(~ Origin, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Private Allele Counts per Origin",
    x = "Allele Count",
    y = "Number of Alleles"
  )

ggsave(paste0("Plots/PA/", m, "_Origin_distribution_of_allele_counts.pdf"), width = 10, height = 8)
ggsave(paste0("Plots/PA/", m, "_Origin_distribution_of_allele_counts.png"), width = 10, height = 8)

#### Collect Fixed Alleles (like filtered_matrix) ####

# Make diploid counts df
diploid_counts <- c(
  Kruger = 16,
  Etosha = 24,
  Kgalagadi = 12
)

origins_filtered_matrix <- origins_private_only_matrix[
  origins_private_only_matrix$Kruger == diploid_counts["Kruger"] |
    origins_private_only_matrix$Etosha == diploid_counts["Etosha"] |
    origins_private_only_matrix$Kgalagadi == diploid_counts["Kgalagadi"],
]

# View the filtered matrix
head(origins_filtered_matrix)
summary(origins_filtered_matrix)
sum(colSums(origins_filtered_matrix != 0)) # sanity check
nrow(origins_filtered_matrix) # should equal above

## Extra sanity check...
# Assume filtered_matrix is your matrix of SNPs x populations
# Convert to data frame
origins_filtered_df <- as.data.frame(origins_filtered_matrix)

# Create an empty list to store results per population
allele_count_freq <- list()

for(pop in colnames(origins_filtered_df)) {
  # Count how many SNPs have each allele count
  freq_table <- origins_filtered_df %>%
    count(Allele_Count = !!sym(pop)) %>%   # count occurrences of each allele count
    arrange(Allele_Count)
  
  # Add a column for the population name
  freq_table$Population <- pop
  
  allele_count_freq[[pop]] <- freq_table
}

# Combine all populations into one data frame
allele_count_freq_df <- bind_rows(allele_count_freq)

# Optional: remove zeros if you don't want them
allele_count_freq_df_no0 <- allele_count_freq_df %>% filter(Allele_Count > 0)

# View
print(allele_count_freq_df_no0) # <- all of these should equal the population_diploid_counts

###############################
# Plot
###############################

library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------------
# 1️⃣ SOURCE: transposed_matrix
# -------------------------------

source_fixed_private_alleles_count <- colSums(filtered_matrix != 0)
source_total_private_alleles_count <- colSums(transposed_matrix != 0)

source_allele_summary <- data.frame(
  Population = colnames(transposed_matrix),
  Fixed_Private_Count = source_fixed_private_alleles_count,
  Total_Private_Count = source_total_private_alleles_count
)

# Add non-fixed private alleles (Total - Fixed)
source_allele_summary$Non_Fixed_Count <- source_allele_summary$Total_Private_Count - source_allele_summary$Fixed_Private_Count

# If there are populations with ZERO fixed alleles
source_allele_summary <- source_allele_summary %>%
  mutate(
    Fixed_Private_Count = as.numeric(ifelse(is.na(Fixed_Private_Count), 0, Fixed_Private_Count)),
    Total_Private_Count = as.numeric(ifelse(is.na(Total_Private_Count), 0, Total_Private_Count))
  )

source_allele_summary

write.table(source_allele_summary, file = paste0(m, "_source_allele_summary.txt"), sep = "\t", row.names = FALSE)

# -------------------------------
# 2️⃣ ORIGIN: origins_private_only_matrix
# -------------------------------

origin_fixed_private_alleles_count <- colSums(origins_filtered_matrix != 0)
origin_total_private_alleles_count <- colSums(origins_private_only_matrix != 0)

origin_allele_summary <- data.frame(
  Population = colnames(origins_private_only_matrix),
  Fixed_Private_Count = origin_fixed_private_alleles_count,
  Total_Private_Count = origin_total_private_alleles_count
)

# Add non-fixed private alleles (Total - Fixed)
origin_allele_summary$Non_Fixed_Count <- origin_allele_summary$Total_Private_Count - origin_allele_summary$Fixed_Private_Count

# If there are populations with ZERO fixed alleles
origin_allele_summary <- origin_allele_summary %>%
  mutate(
    Fixed_Private_Count = as.numeric(ifelse(is.na(Fixed_Private_Count), 0, Fixed_Private_Count)),
    Total_Private_Count = as.numeric(ifelse(is.na(Total_Private_Count), 0, Total_Private_Count))
  )

origin_allele_summary

write.table(origin_allele_summary, file = paste0(m, "_origin_allele_summary.txt"), sep = "\t", row.names = FALSE)

# -------------------------------
# 3️⃣ Combine SOURCE + ORIGIN
# -------------------------------

# Add a "Group" column to distinguish source vs origin
source_allele_summary <- source_allele_summary %>%
  mutate(Group = "Source")

origin_allele_summary <- origin_allele_summary %>%
  mutate(Group = "Origin")

# Combine both
combined_allele_summary <- bind_rows(source_allele_summary, origin_allele_summary)
combined_allele_summary

# Pivot to long format for stacked plotting
combined_allele_long <- combined_allele_summary %>%
  pivot_longer(
    cols = c(Fixed_Private_Count, Non_Fixed_Count),
    names_to = "Condition",
    values_to = "Count"
  )

# Clean up Condition labels for plotting
combined_allele_long$Condition <- recode(
  combined_allele_long$Condition,
  Fixed_Private_Count = "Fixed",
  Non_Fixed_Count = "Non-Fixed"
)

# Preview
combined_allele_long

# -------------------------------
# 4️⃣ Plot
# -------------------------------

# Filter out populations not present in each group
combined_allele_long <- combined_allele_long %>%
  group_by(Group) %>%
  filter(Count > 0) %>%
  ungroup()

ggplot(combined_allele_long, aes(fill = Condition, y = Count, x = Population)) + 
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~Group, scales = "free_x") +   # allow each facet its own populations
  labs(
    #title = "Private Alleles by Population",
    x = "Population",
    y = "Number of Private Alleles"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("Plots/PA/", m, "_Private_Alleles_Origin_and_Source.pdf"), width = 10, height = 10)
ggsave(paste0("Plots/PA/", m, "_Private_Alleles_Origin_and_Source.png"), width = 10, height = 10)

# Filter to include only "Fixed" counts
combined_fixed_counts <- combined_allele_long %>%
  filter(Condition == "Fixed")

# Plot the data
ggplot(combined_fixed_counts, aes(fill = Condition, y = Count, x = Population)) + 
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~Group, scales = "free_x") +   # allow each facet its own populations
  labs(
    #title = "Private Alleles by Population",
    x = "Population",
    y = "Number of Private Alleles"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("Plots/PA/", m, "_Fixed_Alleles_Origin_and_Source.pdf"), width = 10, height = 8)
ggsave(paste0("Plots/PA/", m, "_Fixed_Alleles_Origin_and_Source.png"), width = 10, height = 8)
