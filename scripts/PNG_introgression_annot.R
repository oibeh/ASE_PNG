####################################################
## Description: Genomic annotation of archaic SNPs
####################################################


# =================
# Load libraries
# =================

library(dplyr)
library(tidyverse)
library(valr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(plyranges)



# ===============================================================================
# Calculate odds ratios for the relative enrichment of aSNPs in genomic regions
# ===============================================================================


annotations <- c(
  "hg38_genes_promoters",
  "hg38_genes_exons",
  "hg38_genes_5UTRs",
  "hg38_genes_3UTRs",
  "hg38_genes_introns",
  "hg38_enhancers_fantom"
)

annotatr_annotations <- build_annotations(genome = "hg38", annotations = annotations)

annotated_snps <- annotate_regions(
  regions = aSNPs_haplotypes_v1_pap_fromfilt_gr, #GRanges object of SNPs
  annotations = annotatr_annotations,
  ignore.strand = TRUE
)
archaic_df <- as.data.frame(annotated_snps)


#non-archaic
non_archaic_gr <- GRanges(
  seqnames = non_archaic_df$seqnames,
  ranges = IRanges(start = non_archaic_df$pos, end = non_archaic_df$pos),
  strand = "*"
)

non_archaic_annotated <- annotate_regions(
  regions = non_archaic_gr,
  annotations = annotatr_annotations,
  ignore.strand = TRUE
)
non_archaic_df <- as.data.frame(non_archaic_annotated)

regions <- annotations


#count archaic variants
archaic_counts <- archaic_df %>%
  group_by(annot.type) %>%
  summarise(archaic_count = n()) %>%
  ungroup()

#count non-archaic
genome_counts <- non_archaic_df %>%
  group_by(annot.type) %>%
  summarise(genome_count = n()) %>%
  ungroup()

#combine counts and calculate totals
all_counts <- full_join(archaic_counts, genome_counts, by = "annot.type") %>%
  mutate(
    archaic_count = replace_na(archaic_count, 0),
    genome_count = replace_na(genome_count, 0),
    remainder_archaic = sum(archaic_count) - archaic_count,
    remainder_genome = sum(genome_count) - genome_count
  )

#function to compute bootstrap confidence intervals for odds ratios
bootstrap_or <- function(row, n = 1000) {
  boot_or <- replicate(n, {
    # Resample counts
    archaic_sample <- rbinom(1, size = row$archaic_count + row$remainder_archaic, prob = row$archaic_count / sum(row$archaic_count, row$remainder_archaic))
    genome_sample <- rbinom(1, size = row$genome_count + row$remainder_genome, prob = row$genome_count / sum(row$genome_count, row$remainder_genome))
    # Compute odds ratio
    or <- (archaic_sample / (row$archaic_count + row$remainder_archaic - archaic_sample)) /
          (genome_sample / (row$genome_count + row$remainder_genome - genome_sample))
    or
  })
  c(
    lower = quantile(boot_or, 0.025, na.rm = TRUE),
    upper = quantile(boot_or, 0.975, na.rm = TRUE),
    mean_or = mean(boot_or, na.rm = TRUE)
  )
}

#apply bootstrap to calculate odds ratios and confidence intervals
all_counts <- all_counts %>%
  rowwise() %>%
  mutate(
    or = (archaic_count / remainder_archaic) / (genome_count / remainder_genome),
    ci_lower = bootstrap_or(cur_data())["lower"],
    ci_upper = bootstrap_or(cur_data())["upper"],
    mean_or = bootstrap_or(cur_data())["mean_or"]
  ) %>%
  ungroup()


#perform Fisher's exact test for each region to get p-values
all_counts <- all_counts %>%
  rowwise() %>%
  mutate(
    fisher_p = fisher.test(
      matrix(
        c(archaic_count, remainder_archaic, genome_count, remainder_genome),
        nrow = 2
      )
    )$p.value
  ) %>%
  ungroup()

#annotate significance levels
all_counts <- all_counts %>%
  mutate(
    significance = case_when(
      fisher_p < 0.001 ~ "***",
      fisher_p < 0.01 ~ "**",
      fisher_p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


#plot enrichment/depletion
pdf("genomic_enrichment_oddsratio.pdf", height=5)
ggplot(all_counts, aes(x = reorder(annot.type, -or), y = or)) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "black") +
  geom_text(aes(label = significance), vjust = -0.5, size = 5, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Genomic Region",
    y = "Odds Ratio (Enrichment or Depletion)"
  )
dev.off()
