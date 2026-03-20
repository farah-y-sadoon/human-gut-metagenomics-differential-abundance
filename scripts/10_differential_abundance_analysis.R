"Diversity and differential abundance analysis of human gut microbiome between individuals with a vegan diet (n = 3) and omnivore diet (n = 3)"
# This analysis is based on Lecture 12 and 13 of BINF 6110 - Genomic Methods


# Install and load libraries ----
# if (!require("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# devtools::install_github("jsilve24/ALDEx3")
# install.packages("tidyverse")
# install.packages("phyloseq")
# install.packages("biomformat")
# install.packages("vegan")
# install.packages("viridis")
# BiocManager::install("microbiome")
# install.packages("ggrepel")
# install.packages("patchwork")

library(phyloseq)
library(biomformat)
library(ALDEx3)
library(tidyverse)
library(vegan)
library(viridis)
library(microbiome)
library(ggrepel)
library(patchwork)

# Load data ----
# Set working directory
# Get the name of the current folder
current_dir <- basename(getwd())

# Only change directory if we aren't already in the 'scripts' folder
if (current_dir != "scripts") {
  setwd("./scripts/")
}

# Load metadata
df_metadata <- read_csv("../data/metadata/metadata.csv")

# Read the BIOM file
biom_data <- read_biom("../results/bracken/combined.biom")
physeq <- import_biom(biom_data)

# Create OTU table
otu_table <- as.data.frame(t(otu_table(physeq)))
rownames(otu_table) <- gsub("_report_bracken_species", "", rownames(otu_table))

# Diversity Analysis ----
# Rarefaction ----
png("../figs/03_rarefaction_curve.png", width = 3000, height = 2000, res = 300)
rarecurve <- rarecurve(otu_table, step = 100, label = TRUE, col = 1:6, lty = 1, lwd = 2,
                        xlab = "Sequencing Depth",
                        ylab = "Number of Species")
dev.off()

## Relative Abundance ----
# Calculate relative abundance for each sample
physeq_rel <- transform_sample_counts(physeq,function(x) x / sum(x))
physeq_phy <- tax_glom(physeq_rel, taxrank = "Rank2")

# Create data frame from phyloseq object for plotting
df_phy <- psmelt(physeq_phy)
df_phy$Sample <- gsub("_report_bracken_species", "", df_phy$Sample)
df_phy$Rank2 <- gsub("p__", "", df_phy$Rank2)

# Combine phylum data frame 
df_phy <- df_phy %>% 
  left_join(df_metadata, by = c(Sample = "srr")) %>% 
  select(-c(Rank1, Id)) %>% 
  rename(srr = Sample, otu = OTU, abundance = Abundance, rank2 = Rank2) %>% 
  mutate(rank2 = replace(rank2, rank2 == "", "Unknown"))

# Create relative abundance plot
p4 <- ggplot(df_phy, aes(x = srr, y = abundance, fill = rank2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis(discrete = TRUE, option = "D") + 
  labs(y = "Relative Abundance", x = "Sample (SRR)", fill = "Phylum") + 
  facet_wrap(~diet, scales = "free_x")
ggsave("../figs/04_phylum_abundance.png", plot = p4, width = 10, height = 6.67, dpi = 300, units = "in")
                                      
## Alpha Diversity Analysis ----
### Dominance / Evennness - Berger-Parker ----
# Calculate Berger-Parker Index 
berger_parker <- dominance(physeq, index = "DBP")

# Convert to dataframe and clean up
berger_parker <- berger_parker %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  mutate("Sample (SRR)" = gsub("_report_bracken_species", "", Sample)) %>%
  rename("Berger Parker" = dbp) %>% 
  select("Sample (SRR)", "Berger Parker") %>% 
  remove_rownames()

# Plot Berger Parker (boxplot)
# Prepare data frame for plotting
berger_parker <- berger_parker %>% 
  left_join(df_metadata, by = c("Sample (SRR)" = "srr"))

# Calculate mean to plot
mean_bp <- berger_parker %>%
  group_by(diet) %>%
  summarize(mean_bp = mean(`Berger Parker`))

# Create plot
bp_plot <- ggplot(berger_parker, aes(x = `Sample (SRR)`, y = `Berger Parker`, color = diet, label = `Berger Parker`)) + 
  geom_point(size = 4, alpha = 0.7) +
  scale_color_viridis(discrete = TRUE, option = "D", begin = 0.2, end = 0.7) + 
  geom_hline(data = mean_bp, aes(yintercept = mean_bp, color = diet),
             linetype = "dashed") + 
  geom_text_repel(aes(label = round(`Berger Parker`, 3)),
                  size = 3, max.overlaps = 20) +
  labs(title = "A", x = NULL, y = "Berger-Parker Index", color = "Diet") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../figs/berger_parker_plot.png", plot = bp_plot, width = 10, height = 6.67, dpi = 300, units = "in")

### Information - Shannon ----
# Calculate Shannon Information 
shannon <- estimate_richness(physeq, measures = c("Shannon"))

# Convert to data frame and clean up
shannon <- shannon %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Sample (SRR)") %>% 
  mutate("Sample (SRR)" = gsub("_report_bracken_species", "", `Sample (SRR)`)) %>% 
  remove_rownames()

# Prepare data frame for plotting
shannon <- shannon %>% 
  left_join(df_metadata, by = c("Sample (SRR)" = "srr"))

# Calculate mean to plot
mean_shannon <- shannon %>%
  group_by(diet) %>%
  summarize(mean_shannon = mean(Shannon))

# Plot Shannon Information 
shannon_plot <- ggplot(shannon, aes(x = `Sample (SRR)`, y = Shannon, color = diet, label = Shannon)) + 
  geom_point(size = 4, alpha = 0.7) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  geom_text_repel(aes(label = round(Shannon, 3)), size = 3, max.overlaps = 20) +
  geom_hline(data = mean_shannon, aes(yintercept = mean_shannon, color = diet),
             linetype = "dashed") + 
  labs(title = "B", x = NULL, y = "Shannon Index", color = "Diet") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave("../figs/shannon_plot.png", plot = shannon_plot, width = 10, height = 6.67, dpi = 300, units = "in")

# Save alpha diversity plot and measures
p5 <- bp_plot / shannon_plot + plot_layout(guides = "collect")
ggsave("../figs/05_alpha_diversity_plots.png", plot = p5, width = 10, height = 6.67, dpi = 300, units = "in")

alpha_diversity <- shannon %>%
  left_join(berger_parker, by = c("Sample (SRR)", "diet")) %>%
  select(`Sample (SRR)`, diet, Shannon, `Berger Parker`)

write.csv(alpha_diversity, "../results/R/alpha_diversity.csv", row.names = FALSE)


## Beta Diversity Analysis ----
# Combine metadata with phyloseq object for beta diversity analysis
physeq_beta <- physeq
sample_names(physeq_beta) <- gsub("_report_bracken_species", "", sample_names(physeq_beta))
sample_data(physeq_beta) <- df_metadata %>%
  column_to_rownames("srr")

### PCoA Bray-Curtis ----
ord.pcoa.bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")

bray_pcoa_plot <- plot_ordination(physeq_beta, ord.pcoa.bray, color = "diet") + 
  geom_point(size = 5) + 
  scale_color_viridis_d(option = "D", begin = 0.15, end = 0.65) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.65) +
  labs(title = "A. Bray-Curtis PCoA", color = "Diet") +
  theme_minimal()

### PCoA Jaccard ----
ord.pcoa.jaccard <- ordinate(physeq_beta, method = "PCoA", distance = "jaccard")
jaccard_pcoa_plot <- plot_ordination(physeq_beta, ord.pcoa.jaccard, color = "diet") + 
  geom_point(size = 5) + 
  scale_color_viridis_d(option = "D", begin = 0.15, end = 0.65) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.65) +
  labs(title = "B. Jaccard PCoA") +
  theme_minimal() + 
  theme(legend.position = "none")

# Save beta-diversity plots
p6 <- bray_pcoa_plot + jaccard_pcoa_plot + plot_layout(guides = "collect")
ggsave("../figs/06_beta_diversity_plots.png", plot = p6, width = 10, height = 6.67, dpi = 300, units = "in")

### PERMANOVA ----
# Setup metadata 
metadata_beta <- as(sample_data(physeq_beta), "data.frame")

#### Bray-Curtis PERMANOVA ----
set.seed(123)
permanova_bray <- adonis2(phyloseq::distance(physeq_beta, method = "bray") ~ diet, data = metadata_beta)

# Extract results
bray_result <- as.data.frame(permanova_bray) %>%
  rownames_to_column(var = "Term") %>%
  mutate(Distance = "Bray-Curtis")

#### Jaccard PERMANOVA ----
set.seed(123)
permanova_jaccard <- adonis2(phyloseq::distance(physeq_beta, method = "jaccard") ~ diet, data = metadata_beta)

# Extract results
jaccard_result <- as.data.frame(permanova_jaccard) %>%
  rownames_to_column(var = "Term") %>%
  mutate(Distance = "Jaccard")

# Combine PERMANOVA results and reorder columns
permanova_results <- bind_rows(bray_result, jaccard_result) %>%
  select(Distance, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)

# Save
write.csv(permanova_results, "../results/R/permanova_results.csv", row.names = FALSE)

# Differential Abundance Analysis (DAA) ----
# Prepare matrix as expected by ALDEx3
count_matrix <- as.matrix(otu_table(physeq_beta))
if(!taxa_are_rows(physeq_beta)) {
  count_matrix <- t(count_matrix) # transpose table if taxa are not rows
  } 

set.seed(123)
aldex_out <- aldex(Y = count_matrix, 
                   X = ~diet, 
                   data = metadata_beta,
                   method = "lm",
                   scale = clr.sm,
                   p.adjust.method = "BH",
                   test = "t.HC3")

# Pull out relevant results
aldex_out_summary <- summary(aldex_out)

# Join with phyloseq object to get taxa names for plotting
df_aldex_species <- as.data.frame(tax_table(physeq_beta)) %>%
  rownames_to_column("entity") %>%
  mutate(Rank6 = gsub("g__", "", Rank6),
         Rank7 = gsub("s__", "", Rank7),
         species = paste(Rank6, Rank7)) %>% 
  select(entity, species)

aldex_out_summary <- aldex_out_summary %>% 
  left_join(df_aldex_species, by = c(entity = "entity"))

# Plot ALDEx3 differential abundance results 
p7 <- aldex_out_summary %>%
  filter(abs(estimate) > 5) %>%
  ggplot(aes(x = estimate, y = reorder(species, estimate))) +
  geom_point(aes(color = p.val.adj < 0.05), size = 3) +
  geom_errorbar(aes(xmin = estimate - std.error,
                    xmax = estimate + std.error)) +
  geom_vline(xintercept = 0, color = "red") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red"),
                     name = "Adjusted p-value",
                     labels = c("FALSE" = "≥ 0.05", "TRUE" = "< 0.05")) +
  labs(x = "Estimated CLR Abundance Difference (Vegan vs. Omnivore)",
       y = "Species") +
  theme_minimal()

# Save ALDEx3 results and plot
write.csv(aldex_out_summary, "../results/R/aldex3_results.csv", row.names = FALSE)
ggsave("../figs/07_aldex3_daa.png", plot = p7, width = 10, height = 8, dpi = 300, units = "in")
