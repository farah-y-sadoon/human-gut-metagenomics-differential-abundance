# Exploring The Relationship Between Diet and Gut Microbiome Composition Through Shotgun Metagenomics

## Introduction
Metagenomics is a powerful approach for deepening our understanding of the microbial diversity and community dynamics of the human gut (1-3). Metagenomic approaches aim to characterize the taxonomic composition of samples. In the context of the human gut, this primarily means cataloging and pinpointing the function of bacterial species associated with health and disease (1, 2). Within metagenomics, metabarcoding and shotgun metagenomics are two strategies used to understand the composition and functionality of microbial communities. Metabarcoding relies on the targeted sequencing of amplicons to detect taxa with relatively low coverage; however, it does not provide a direct link to functional profiling and is highly reliant on the chosen amplicon for species recovery (4). Shotgun metagenomics sequences all DNA within samples, allowing for the recovery of rare taxa and direct functional profiles, although it requires substantially higher sequencing depth and produces false positive results (4). For both methods, classification and statistical software and database completeness are major limiting factors, meaning that choosing appropriate tools is critical for an effective study (5-7). This analysis focuses on exploring shotgun metagenomics workflows for investigating the taxonomic composition of the human gut microbiome.

Some of the most common k-mer-based tools for taxonomic classification are Kraken2 (8), KrakenUniq (9), and CLARK (10). When choosing a classifier, precision, recall, accuracy, and false positive rate are important metrics to consider (5). In a benchmarking study, Ye and colleagues found that all three classifiers had a median area under the precision-recall curve (AUPR) score of about 0.95 and L2 values below 0.1. These results underscore the tools’ high precision and recall, as well as accuracy against ground-truth data. While all three tools perform well on the described metrics, computational resources must also be considered. Kraken2 and CLARK ran in a similar amount of time on a dataset with approximately 5.7 million reads (~10¹ minutes), while KrakenUniq ran considerably slower (~10² minutes). Additionally, KrakenUniq requires hundreds of gigabytes of memory. While all three tools are similar in their performance for classification, KrakenUniq is the superior choice due to its reduced false positive and misclassification rate; however, given resource constraints, Kraken2 with Bracken (11) is a suitable alternative (5). For these reasons, this study uses Kraken2 with Bracken for classification.

Database constraints represent a major limiting factor in metagenomic analyses (7). Some available databases for Kraken2 are the standard-16, standard, and nt-core databases. At varying confidence scores, each database performs differently in its ability to identify reads at taxonomic levels (12). Notably, the standard-16 database classifies 0% of reads at 40% confidence, while the standard and nt-core databases classify 80% and 95% at the same confidence level, respectively. Interestingly, it seems that confidence score has a larger impact on precision, recall, and F1 score than the database itself (12). Moreover, computational constraints should be considered, as larger databases will contain more records for accuracy in classification, but require more memory and storage. Here, the Kraken2 standard database was used with a confidence score of 0.15 to balance classification rates, precision, and recall.

Selecting a statistical tool for differential abundance analysis is also crucial. ANCOM-BC2 (13) and ALDEx-2 (14) are both recommended; they give consistent results and help control false positives (6). Recent benchmarking studies show that ALDEx2 is more consistent. For example, when results between exploratory and validation datasets were compared (15), ANCOM-BC2 had higher conflict rates (3%) and lower replication rates (35%). ALDEx2 performed better, with only 0.1% conflict and 79% replication at a false discovery rate (FDR) of 0.05. ALDEx2 is considered conservative, but its successor, ALDEx3 (16), builds on the same framework and adds improvements. ALDEx3 was chosen for this workflow to maintain control of false positives.

Several diseases, such as type 2 diabetes, cardiovascular disease, and irritable bowel disease, have been linked to the composition and functionality of the gut microbiome (17). While the complexities of human gut microbiota are still being explored, diet has emerged as a major contributing factor (17). In particular, Westernized, low-fiber diets are linked to reduced microbial diversity and loss of fibre-degrading taxa (17). Studying gut microbial taxa in people with different diets can help inform recommendations to reduce chronic health issues related to the microbiome. Shotgun metagenomic studies across dietary groups provide insights into specific taxa for further study (18).

The study described here explores a bioinformatic workflow for processing shotgun metagenomics data, and applies this workflow to examine differences in gut microbiome composition between individuals with vegan and omnivore diets.

## Methods
### Computational Resources

Computational workflows, including data acquisition, SRA conversion, quality assessment, read trimming and filtering, taxonomic classification, and abundance estimation, were executed on the Digital Research Alliance of Canada’s Fir cluster (19), with most steps submitted as SLURM jobs. Downstream analyses, including diversity metrics and differential abundance testing, were performed locally on a MacBook Pro (M4 architecture). Mamba v2.4.0 (20) managed virtual environments and software dependencies throughout the workflow (see PIPELINE.md). Git was used for version control (21).

### Data Acquisition

The Kraken2 standard reference v20251015 database was obtained from the Kraken2 documentation index zone (22). The `prefetch` and `fasterq-dump` functions from SRA Toolkit v3.0.9 (23) were used to download and convert SRA files to FASTQ format.


### Quality Control

Quality assessment of reads was conducted for each sample with FastQC v0.12.1 (24) and consolidated into a single report with MultiQC v1.33 (25) before filtering reads. Fastp v1.0.1 (26) was run in parallel mode using `parallel.py` with `--detect_adapter_for_pe` to remove remaining adapters, `--length_required 50` to ensure reads were larger than Kraken2 k-mer size, and `qualified_quality_phred 20` with `--unqualified_percent_limit 20` to remove reads with more than 20 percent of bases with Phred scores below 20. Since Fastp provides its own summary report, MultiQC was not run after trimming and filtering.

### Taxonomic Classification

Kraken2 v2.1.6 (8) was run for each sample, with `confidence 0.15` to reduce the false positive rate, and `--paired` to indicate that each sample includes forward and reverse reads. Bracken v3.0 (11) was run to perform abundance re-estimation with `-r 150` to reflect the length of raw reads, and `-l S` to include only species-level estimates. To convert individual Braken reports to a combined BIOM object for downstream diversity analyses, kraken-biom v1.2.0 (27) was used with `--json` to ensure compatibility with R v4.5.1 (28).

### Diversity and Differential Abundance Analysis

The combined BIOM file was loaded into R and converted into a phyloseq object with the phyloseq package v1.52.0 (29). After creating an operational taxonomic unit (OTU) table, rarefaction was conducted with the vegan v2.7-3 `rarecurve` function (30) to assess sampling completeness. Alpha diversity was measured with the Berger-Parker index to assess dominance, and Shannon index to determine richness and evenness. These metrics were chosen to follow recommendations to include alpha diversity metrics measuring multiple aspects (31). The `dominance` function from the microbiome package v1.30.0 was used with `index = “DBP”` to perform Berger-Parker calculations, and the `estimate_richness` function from phyloseq was used with `measures = c("Shannon")` to calculate Shannon index values. Bray-Curtis dissimilarity and Jaccard similarity were calculated with the `ordinate` function from phyloseq, and PERMANOVA was conducted with the `adonis2` function from vegan to determine the statistical significance of species composition between vegan and omnivore samples. Finally, differential abundance analysis was conducted with ALDEx3 v1.0.2 (16) to determine whether bacterial taxa were significantly enriched or depleted between vegan and omnivore samples, accounting for compositionality and variability in sequencing depth. The “lm” method was specified to assess the linear relationship between diet and abundance with the Benjamini-Hochberg correction. The scale method was set to “clr.sm” and the test parameter was set to "t.HC3" as recommended by the ALDEx3 documentation for small sample sizes (16).

### Visualizations

Visualizations were generated using `rarecurve` from vegan, and `ggplot` from ggplot2 v4.0.2 (32). 

### Pipeline
```mermaid
flowchart LR
    B["Data Acquisition (Raw reads: SRA Toolkit, Database: Kraken2 DB)"]
    B --> C["Quality Control (FastQC, MultiQC, Fastp)"]
    C --> D["Taxonomic Classification (Kraken2, Bracken, kraken-biom)"]
    D --> E["Diversity & Differential Abundance (R, phyloseq, vegan, microbiome, ALDEx3)"]
    E --> F["Visualization (ggplot2, vegan)"]
```
Figure 1. Workflow used for data acquisition, taxonomic classification and diversity analyses to assess differences in gut microbiome composition in human gut samples of vegan and omnivore groups.

## Results
### Statistics after read trimming and filtering show minimal loss of information and improved quality

Trimming and filtering with Fastp (26) led to an increased proportion of high-quality Phred scores, where post-filtering Q20 rates ranged from 94.7% to 98.0% and Q30 rates ranged from 85.6% to 94.0% across all samples (Table 1). Additionally, GC content was preserved, indicating that there was no notable bias introduced during the filtration and trimming process. As expected by Illumina reads, quality scores were lower for reverse reads across all samples, indicated by the decreased number of total reads after processing (33, Table 1).

Table 1. **Read quality metrics before and after adapter trimming and quality filtering with fastp v1.0.1.** Per-sample summary statistics for all six samples (forward and reverse reads reported separately), including total reads, total bases, Q20 rate, Q30 rate, and GC content. Filtering was performed in paired-end mode with a minimum read length of 50 bp to satisfy Kraken2 k-mer requirements, and a Phred quality threshold of Q20 with a maximum 20% unqualified bases per read. Quality metrics improved consistently across all samples following filtering.

| Sample Reads | Total Reads (Before) | Total Reads (After) | Total Bases (Before) | Total Bases (After) | Q20 Rate (Before) | Q20 Rate (After) | Q30 Rate (Before) | Q30 Rate (After) | GC Content (Before) | GC Content (After)
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -- |
SRR8146972 (Forward) | 27.31M | 26.40M | 4.10G | 3.95G | 96.90% | 97.74% | 91.85% | 93.34% | 52.40% | 52.38%
SRR8146972 (Reverse) | 27.31M | 24.60M | 4.10G | 3.69G | 93.90% | 96.91% | 86.43% | 91.31% | 52.48% | 52.35%
SRR8146973 (Forward) | 34.58M | 33.79M | 5.19G | 5.06G | 97.42% | 97.99% | 92.94% | 93.99% | 55.36% | 55.36%
SRR8146973 (Reverse) | 34.58M | 31.83M | 5.19G | 4.77G | 94.73% | 97.12% | 88.12% | 92.05% | 55.46% | 55.36%
SRR8146974 (Forward) | 35.42M | 34.59M | 5.31G | 5.18G | 97.39% | 97.98% | 92.85% | 93.93% | 56.26% | 56.27%
SRR8146974 (Reverse) | 35.42M | 32.56M | 5.31G | 4.87G | 94.66% | 97.07% | 87.95% | 91.92% | 56.36% | 56.27%
SRR8146975 (Forward) | 35.60M | 34.54M | 5.34G | 5.17G | 97.06% | 97.82% | 92.27% | 93.62% | 47.23% | 47.20%
SRR8146975 (Reverse) | 35.60M | 32.45M | 5.34G | 4.86G | 94.35% | 97.08% | 87.36% | 91.80% | 47.32% | 47.15%
SRR8146976 (Forward) | 28.58M | 25.62M | 4.29G | 3.84G | 92.88% | 95.42% | 83.25% | 87.50% | 52.71% | 52.46%
SRR8146976 (Reverse) | 28.58M | 23.23M | 4.29G | 3.48G | 89.23% | 94.66% | 77.18% | 85.58% | 52.82% | 52.34%
SRR8146977 (Forward) | 39.24M | 37.99M | 5.89G | 5.68G | 97.00% | 97.81% | 92.14% | 93.60% | 47.93% | 47.90%
SRR8146977 (Reverse) | 39.24M | 35.70M | 5.89G | 5.36G | 94.35% | 97.11% | 87.37% | 91.89% | 48.01% | 47.80%

### Rarefaction curves indicate sampling completeness
Since each sample was sequenced at different depths, rarefaction curves were generated to assess sampling completeness. Variation in species richness between samples ranges from a minimum of approximately 290 species to a maximum of approximately 375 species (Figure 2). All samples approach an asymptote well before maximum depth, so any differences in composition are unlikely to be caused by sampling bias (Figure 2).

![Figure 2](./figs/02_rarefaction_curve.png)
Figure 2. **Rarefaction curves indicate adequate sequencing depth across all samples.** Species richness (number of detected species) plotted as a function of sequencing depth for each of the six samples (SRR8146972–SRR8146977). All curves approach an asymptote well before maximum sequencing depth (~2.5 × 10⁷ reads). Sample SRR8146975 exhibits the lowest species richness (~290 species), while SRR8146973 shows the highest (~375 species).

### Vegan and omnivore diet-groups display differences in phylum-level relative abundance 

Investigation of phylum-level relative abundance revealed that all samples show a high proportion of Bacillota and Bacteroidota (Figure 3). Additionally, phylum Actinomycetota comprises a notable fraction of the relative abundance for samples SRR8146976, SRR8146973, and SRR8146974, with two of the three belonging to the vegan group (Figure 3). The same vegan samples also showed higher abundance of Verrucomicrobiota. In particular, relative abundance for sample SRR8146974 is comprised of approximately 20% Verrucomicrobiota (Figure 3). The remaining phyla are present at small proportions. 

![Figure 3](./figs/03_phylum_abundance.png)
Figure 3. **Phylum-level community composition shows differences between omnivore and vegan diet groups.** Stacked bar plots showing relative abundance at the phylum level for each of the six samples (SRR8146972–SRR8146977). All samples have a high proportion of Bacillota, Bacteroidota. Actinomycetota comprises a notable fraction across some samples, and appears at a higher proportion overall in the vegan group. The vegan group also showed more Verrucomicrobiota abundance, most prominently in sample SRR8146974 which showed approximately 20% relative abundance for this phylum. Minor phyla collectively account for a small proportion of reads across all samples.

### Berger-Parker and Shannon indices indicate differences in species dominance but similar overall diversity between vegan and omnivore groups
Exploration of alpha diversity metrics showed that the omnivore diet group (mean = ~0.22) had higher species dominance than the vegan diet group when measured by the Berger-Parker index (mean = ~0.15, Figure 4A). Conversely, Shannon index values revealed similar overall diversity for both groups (mean = ~3.27, Figure 4B); however, vegan samples displayed a larger range of Shannon index values, indicating higher variability of species-level diversity in this group (Figure 4B).

![Figure 4](./figs/04_alpha_diversity_plots.png)
Figure 4. **Alpha diversity indices reveal higher species dominance in omnivores but similar overall diversity between groups.** Dot plots showing alpha diversity indices for omnivore (purple) and vegan (green) samples, with dashed lines indicating group means. (A) Berger-Parker index, a measure of dominance, was higher in omnivores (mean ~0.22) than vegans (mean ~0.15). (B) Shannon index, a measure for overall diversity accounting for both richness and evenness, was similar between groups (~3.27).

### Bray-Curtis and Jaccard indices reveal high within-group variation and no clear separation by diet group
Following the analysis of alpha diversity metrics, beta diversity analysis was performed to assess differences in community composition. Bray-Curtis distance incorporates both species presence and abundance, while Jaccard distance only considers species presence and absence. Both beta diversity measures were used to calculate dissimilarity between all samples. Bray-Curtis PCoA explains 37% and 28.5% of variance in PC1 and PC2, respectively (Figure 5A), while Jaccard PCoA explains 30.1% and 26.3% of variance in PC1 and PC2, respectively (Figure 5B). Within-group variation appears comparable to or greater than between-group variation, with no clear clustering by diet group observed by either metric (Figure 5). PERMANOVA further confirmed no significant difference in community composition between diet groups for either distance metric (Table 2). 

![Figure 5](./figs/05_beta_diversity_plots.png)
Figure 5. **Beta diversity metrics show weak separation between omnivore and vegan diet groups.** Principal Coordinates Analysis (PCoA) plot showing community composition dissimilarity within and between omnivore (purple) and vegan (green) diet groups, using two distance metrics. (A) Bray-Curtis PCoA, which accounts for species presence and abundance, explains 37% and 28.5% of variance in PC1 and PC2, respectively. (B) Jaccard PCoA, which only considers species presence or absence, explains 30.1% and 26.3% of variance in PC1 and PC2, respectively. Both metrics reveal high within-group variation and no clear clustering by diet group.

Table 2. **PERMANOVA results for Bray-Curtis and Jaccard distance matrices.** PERMANOVA was performed using the `adonis2` function in the vegan package to test for significant differences in community composition between diet groups.
| Distance Metric | df | R² | F | p-value |
| --- | --- | --- | --- | --- |
| Bray-Curtis | 1 | 0.246 | 1.303 | 0.300 |
| Jaccard | 1 | 0.234 | 1.222 | 0.200 |

### Differential abundance analysis of community composition reveals weak trends in estimated CLR between groups

ALDEx3 was used to calculate Centered Log-Ratio (CLR) estimates for species abundance counts across all samples. The CLR transformation is a normalization method that is applied to compositional data, allowing for the comparison of relative abundance between samples (34). After transformation, there were no significant differences in species abundance between the vegan and omnivore groups; however, some directional trends still emerge (Figure 6). Most notably, some *Prevotella* species show higher relative abundance in the omnivore group, while *Akkermansia muciniphila* has higher relative abundance in the vegan group (Figure 6).

![Figure 6](./figs/06_aldex3_daa.png)
Figure 6. **Estimated CLR abundance values show directional trends but no statistically significant differences in composition between vegan and omnivore diet groups.** Forest plot of estimated CLR (Centered Log-Ratio) abundance differences (vegan vs. omnivore) for the top species by effect size, as determined by ALDEx3 differential abundance analysis. Each point represents a species estimate with 95% confidence intervals. The red vertical line indicates no difference (CLR = 0). Positive values indicate higher relative abundance in vegans, while negative values indicate higher relative abundance in omnivores. All species have adjusted p-values ≥ 0.05 (grey points), indicating no statistically significant differential abundance after Benjamini-Hochberg correction. Notably, most *Prevotella* species trend higher toward omnivores, while *Akkermansia muciniphila* trends higher toward vegans.

## Discussion
This study investigated gut microbiome community composition in individuals with vegan and omnivore diets using a shotgun metagenomics workflow. The pipeline successfully classified over 10 million reads per sample at the species level, and while no statistically significant differences were detected between diet groups, several trends emerged that are consistent with existing literature and warrant further investigation.

The investigation of phylum-level relative abundance between groups showed that Bacillota (also known as Firmicutes) and Bacteroidota were the most represented across all groups (Figure 3). This is consistent with our current understanding that these two phyla represent the highest proportion of relative abundance in the human gut microbiome (35). Actinomycetota (also known as Actinobacteria), which is considered a relevant minority phylum (36, 37), was also noticeably abundant in samples, albeit mostly in the vegan group (Figure 3). Interestingly, there have been contrasting results regarding the effect of diet on the presence of Actinmycetota, with some studies reporting increases among vegan samples, and others reporting a decrease (38). This result suggests that factors other than diet may influence the presence of this phylum. The increased representation of Verrucomicrobiota in vegans is supported by other work investigating the effects of diet on microbiome composition (39). Most research on Verrucomicrobiota investigate its ecological and environmental roles, but some work suggests that *Akkermansia muciniphila* can help to protect against chronic illnesses such as inflammatory bowel disease and type two diabetes (40).

The exploration of alpha diversity revealed that dominance, measured by the Berger-Parker index, was higher in the omnivore group. While the difference was not confirmed with statistical testing, this result suggests that species composition is dominated by one or a few species for the omnivore group, and more evenly spread for the vegan group (Figure 4). No direct link between dominance and diet was identified in the existing literature; however, differences in alpha diversity between vegan and omnivore diet groups have been linked to dietary intake (41). Future work may benefit from explicitly investigating differences in dominance for diet groups, as it may provide insights into the effects of dietary intake on microbiome composition.

Bray-Curtis and Jaccard distances were used to assess between-sample differences in community composition. PERMANOVA analysis provided no evidence of a significant difference between the vegan and omnivore groups (Figure 5, Table 2). A notable limitation of this study was the small sample size for each group (n = 3), as the statistical power of such a small study requires substantial differences between groups to reach significance. Although no significant separation between dietary groups was detected, this may reflect insufficient statistical power rather than a true absence of differences, as significant differences in gut microbiome composition have been linked to diet (41).

To identify specific taxa contributing to compositional differences between groups, differential abundance analysis was conducted using ALDEx3. The results indicated no significant differences between vegan and omnivore groups for CLR estimates of species abundance (Figure 6). Nevertheless, directional trends demonstrated potential differences between groups, which may reflect the low statistical power of the small sample size used here (Figure 6). Specifically, most *Prevotella* species showed higher relative abundance in the omnivore group. Similarly to Fackelman and colleagues’ recent study, *Prevotella copri* was not a strong signature of vegan diets (41). The increased relative abundance of *Prevotella* species in the omnivore group observed here warrants further investigation in a larger cohort to determine whether this reflects a consistent dietary indicator or an artifact of the small sample size. Furthermore, *Akkermansia muciniphila*, a species within Verrucomicrobiota, showed higher relative abundance in the vegan group (Figure 6). *A. muciniphila* is considered a promising probiotic candidate, as it has been linked to improved gut health and an inverse relationship with chronic conditions such as inflammatory bowel disease and type 2 diabetes (42). Future functional experiments would help elucidate the biological mechanisms underlying its abundance and whether diet is an important factor for its colonization and persistence in the human gut.

Overall, this study applied a shotgun metagenomics workflow to investigate gut microbiome composition in vegan and omnivore diet groups. While no statistically significant differences were detected across any of the diversity or differential abundance analyses, consistent trends emerged across multiple metrics. Elevated Verrucomicrobiota and *Akkermansia muciniphila* in vegans, alongside higher species dominance and enrichment of *Prevotella* species in omnivores, suggest that diet may influence gut microbiome composition in ways that are biologically meaningful but require larger sample sizes to detect statistically. Future studies with greater statistical power and functional profiling would help to clarify the relationship between dietary patterns and gut microbiome composition, and further explore the potential health implications of diet-associated taxa.

## References
1. Wang, W.-L., Xu, S.-Y., Ren, Z.-G., Tao, L., Jiang, J.-W., & Zheng, S.-S. (2015). Application of metagenomics in the human gut microbiome. World Journal of Gastroenterology, 21(3), 803. https://doi.org/10.3748/wjg.v21.i3.803

2. Qin, J., Li, R., Raes, J., Arumugam, M., Burgdorf, K. S., Manichanh, C., Nielsen, T., Pons, N., Levenez, F., Yamada, T., Mende, D. R., Li, J., Xu, J., Li, S., Li, D., Cao, J., Wang, B., Liang, H., Zheng, H., … Wang, J. (2010). A human gut microbial gene catalogue established by metagenomic sequencing. Nature, 464(7285), 59–65. https://doi.org/10.1038/nature08821 

3. Gill, S. R., Pop, M., DeBoy, R. T., Eckburg, P. B., Turnbaugh, P. J., Samuel, B. S., Gordon, J. I., Relman, D. A., Fraser-Liggett, C. M., & Nelson, K. E. (2006). Metagenomic analysis of the human distal gut microbiome. Science, 312(5778), 1355–1359. https://doi.org/10.1126/science.1124234 

4. Durazzi, F., Sala, C., Castellani, G., Manfreda, G., Remondini, D., & De Cesare, A. (2021). Comparison between 16S rrna and shotgun sequencing data for the taxonomic characterization of the gut microbiota. Scientific Reports, 11(1). https://doi.org/10.1038/s41598-021-82726-y 

5. Ye, S. H., Siddle, K. J., Park, D. J., & Sabeti, P. C. (2019). Benchmarking metagenomics tools for taxonomic classification. Cell, 178(4), 779–794. https://doi.org/10.1016/j.cell.2019.07.010 

6. Nearing, J. T., Douglas, G. M., Hayes, M. G., MacDonald, J., Desai, D. K., Allward, N., Jones, C. M., Wright, R. J., Dhanani, A. S., Comeau, A. M., & Langille, M. G. (2022). Microbiome differential abundance methods produce different results across 38 datasets. Nature Communications, 13(1). https://doi.org/10.1038/s41467-022-28034-z 

7. Chorlton, S. D. (2024). Ten common issues with reference sequence databases and how to mitigate them. Frontiers in Bioinformatics, 4. https://doi.org/10.3389/fbinf.2024.1278228   

8. Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1). https://doi.org/10.1186/s13059-019-1891-0 

9. Breitwieser, F. P., Baker, D. N., & Salzberg, S. L. (2018). Krakenuniq: Confident and fast metagenomics classification using unique K-Mer counts. Genome Biology, 19(1). https://doi.org/10.1186/s13059-018-1568-0 

10. Ounit, R., Wanamaker, S., Close, T. J., & Lonardi, S. (2015). CLARK: Fast and accurate classification of metagenomic and genomic sequences using discriminative K-MERS. BMC Genomics, 16(1). https://doi.org/10.1186/s12864-015-1419-2

11. Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: Estimating species abundance in metagenomics data. PeerJ Computer Science, 3. https://doi.org/10.7717/peerj-cs.104 

12. Liu, Y., Ghaffari, M. H., Ma, T., & Tu, Y. (2024). Impact of database choice and confidence score on the performance of taxonomic classification using kraken2. aBIOTECH, 5(4), 465–475. https://doi.org/10.1007/s42994-024-00178-0 

13. Lin, H., & Peddada, S. D. (2020). Analysis of compositions of microbiomes with bias correction. Nature Communications, 11(1). https://doi.org/10.1038/s41467-020-17041-7 

14. Fernandes, A. D., Reid, J. N., Macklaim, J. M., McMurrough, T. A., Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of high-throughput sequencing datasets: Characterizing RNA-seq, 16S rrna gene sequencing and selective growth experiments by compositional data analysis. Microbiome, 2(1). https://doi.org/10.1186/2049-2618-2-15 

15. Pelto, J., Auranen, K., Kujala, J. V., & Lahti, L. (2025). Elementary methods provide more replicable results in microbial differential abundance analysis. Briefings in Bioinformatics, 26(2). https://doi.org/10.1093/bib/bbaf130 

16. Justin Silverman (2026). ALDEx3: Linear Models for Sequence Count Data. R package version 1.0.1, https://cran.r-project.org/web/packages/ALDEx3 

17. De Filippis, F., Vitaglione, P., Cuomo, R., Berni Canani, R., & Ercolini, D. (2018). Dietary interventions to modulate the gut microbiome—how far away are we from Precision Medicine. Inflammatory Bowel Diseases, 24(10), 2142–2154. https://doi.org/10.1093/ibd/izy080

18. De Filippis, F., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct genetic and functional traits of human intestinal prevotella copri strains are associated with different habitual diets. Cell Host &amp; Microbe, 25(3). https://doi.org/10.1016/j.chom.2019.01.004 

19. Government of Canada. Digital Research Alliance of Canada. (2023). https://alliancecan.ca  

20. Mamba-org. (n.d.). Mamba: The fast cross-platform package manager [Computer software]. GitHub. https://github.com/mamba-org/mamba 

21. Sadoon, F. (2026). human-gut-metagenomics-differential-abundance [Source code]. GitHub. https://github.com/farah-y-sadoon/human-gut-metagenomics-differential-abundance 

22. Langmead, B. (n.d.). Index zone. Index zone by BenLangmead. https://benlangmead.github.io/aws-indexes/k2 

23. U.S. National Library of Medicine. (n.d.). 01. downloading SRA Toolkit · NCBI/SRA-Tools Wiki · github. National Center for Biotechnology Information. https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software 

24. Babraham Bioinformatics. (n.d.). FastQC. FastQC a quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 

25. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354 

26. Chen, S. (2025). Fastp 1.0: An ultra‐fast all‐round tool for FASTQ data quality control and preprocessing. iMeta, 4(5). https://doi.org/10.1002/imt2.70078 

27. Dabdoub, SM (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software]. Available at https://github.com/smdabdoub/kraken-biom

28. R Core Team. (n.d.). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing. R. https://www.r-project.org/ 

29. McMurdie, P. J., & Holmes, S. (2013). Phyloseq: An R package for reproducible interactive analysis and graphics of Microbiome Census Data. PLOS ONE. https://doi.org/10.1371/journal.pone.0061217 

30. Oksanen, J., et al. (2026). vegan: Community Ecology Package. R package version 2.7-3. https://vegandevs.github.io/vegan/ 

31. Cassol, I., Ibañez, M., & Bustamante, J. P. (2025, January 3). Key features and guidelines for the application of Microbial Alpha Diversity Metrics. Nature News. https://www.nature.com/articles/s41598-024-77864-y 

32. Wickham, H. (n.d.). GGPLOT2. SpringerLink. https://link.springer.com/book/10.1007/978-3-319-24277-4 

33. Kwon, S., Park, S., Lee, B., & Yoon, S. (2013). In-depth analysis of interrelation between quality scores and real errors in Illumina reads. Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE Engineering in Medicine and Biology Society. Annual International Conference. https://pubmed.ncbi.nlm.nih.gov/24109767/

34. Quinn, T. P., Erb, I., Gloor, G., Notredame, C., Richardson, M. F., & Crowley, T. M. (2019). A field guide for the compositional analysis of any-omics data. GigaScience. https://pmc.ncbi.nlm.nih.gov/articles/PMC6755255/ 

35. Qin, J., Li, R., Raes, J., Arumugam, M., Burgdorf, K. S., Manichanh, C., Nielsen, T., Pons, N., Levenez, F., Yamada, T., Mende, D. R., Li, J., Xu, J., Li, S., Li, D., Cao, J., Wang, B., Liang, H., Zheng, H., … Wang, J. (2010). A human gut microbial gene catalogue established by metagenomic sequencing. Nature News. https://www.nature.com/articles/nature08821 

36. Segata, N., Haake, S. K., Mannon, P., Lemon, K. P., Waldron, L., Gevers, D., Huttenhower, C., & Izard, J. (2012). Composition of the adult digestive tract bacterial microbiome based on seven mouth surfaces, tonsils, throat and stool samples - genome biology. SpringerLink. https://link.springer.com/article/10.1186/gb-2012-13-6-r42 

37. Cecilia, B., Riccardo, L. L., Gianenrico, R., Giulia, G., Vincenzo, C., & Antonio, G. (2018). Actinobacteria: A relevant minority for the maintenance of gut homeostasis: Article information: J-global. Digestive and Liver Disease. https://www.sciencedirect.com/science/article/abs/pii/S159086581830210X 

38. Soldán, M., Argalášová, Ľ., Hadvinová, L., Galileo, B., & Babjaková, J. (2024). The effect of dietary types on gut microbiota composition and development of non-communicable diseases: A narrative review. Nutrients. https://pmc.ncbi.nlm.nih.gov/articles/PMC11434870/ 

39. Losno, E. A., Sieferle, K., Perez-Cueto, F. J. A., & Ritz, C. (2021). Vegan diet and the gut microbiota composition in healthy adults. MDPI. https://www.mdpi.com/2072-6643/13/7/2402 

40. Rodrigues, V. F., Elias-Oliveira, J., Pereira, Í. S., Pereira, J. A., Barbosa, S. C., Machado, M. S. G., & Carlos, D. (2022). Akkermansia muciniphila and Gut Immune System: A Good Friendship That Attenuates Inflammatory Bowel Disease, Obesity, and Diabetes. Frontiers in immunology, 13, 934695. https://doi.org/10.3389/fimmu.2022.934695  

41. Fackelmann, G., Manghi, P., Carlino, N., Heidrich, V., Piccinno, G., Ricci, L., Piperni, E., Arrè, A., Bakker, E., Creedon, A. C., Francis, L., Capdevila Pujol, J., Davies, R., Wolf, J., Bermingham, K. M., Berry, S. E., Spector, T. D., Asnicar, F., & Segata, N. (2025). Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. Nature microbiology. https://pmc.ncbi.nlm.nih.gov/articles/PMC11726441/ 

42. Pellegrino, A., Coppola, G., Santopaolo, F., Gasbarrini, A., & Ponziani, F. R. (2023). Role of akkermansia in human diseases: From causation to therapeutic properties. Nutrients. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10142179/ 
