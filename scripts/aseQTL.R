#################################
## Description: aseQTL analysis
#################################


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
library(rstatix)
library(ggpubr)


# =================================
# aseQTL rank permutation testing
# =================================

#extract the lead snps first (rank 1), then rank 2, etc.
rank_subsets <- list()

for (i in 1:5) {
  rank_subsets[[i]] <- raw_cht2_gr_df2_filt2_padj05_aseqtls %>%
    group_by(SYMBOL) %>%
    arrange(P.ADJ) %>%
    slice(i) %>%
    ungroup() %>%
    as.data.frame()
}

#find intersection btw aseQTL ranks and archaic SNPs (by ancestry)
aSNPs_haplotypes_v1_pap_fromfilt2_DEN <- aSNPs_haplotypes_v1_pap_fromfilt2[aSNPs_haplotypes_v1_pap_fromfilt2$refined_ancestry == "Denisova",]
aSNPs_haplotypes_v1_pap_fromfilt2_NEANDPAP <- aSNPs_haplotypes_v1_pap_fromfilt2[aSNPs_haplotypes_v1_pap_fromfilt2$refined_ancestry == "Neand_PAP_uniq",]
aSNPs_haplotypes_v1_pap_fromfilt2_NEAND1KG <- aSNPs_haplotypes_v1_pap_fromfilt2[aSNPs_haplotypes_v1_pap_fromfilt2$refined_ancestry == "Neand_1KG",]

rank1_subset_overlap1 <- filter(aSNPs_haplotypes_v1_pap_fromfilt2_DEN, SNP %in% unique(rank1_subset$SNP))
rank1_subset_overlap2 <- filter(aSNPs_haplotypes_v1_pap_fromfilt2_NEANDPAP, SNP %in% unique(rank1_subset$SNP))
rank1_subset_overlap3 <- filter(aSNPs_haplotypes_v1_pap_fromfilt2_NEAND1KG, SNP %in% unique(rank1_subset$SNP))

#create lists to store the filtered results
rank_subset_overlap1 <- list()
rank_subset_overlap2 <- list()
rank_subset_overlap3 <- list()

#loop through ranks 1-5
for (i in 1:5) {
  rank_subset_overlap1[[i]] <- filter(aSNPs_haplotypes_v1_pap_fromfilt2_DEN, SNP %in% unique(rank_subsets[[i]]$SNP))
  rank_subset_overlap2[[i]] <- filter(aSNPs_haplotypes_v1_pap_fromfilt2_NEANDPAP, SNP %in% unique(rank_subsets[[i]]$SNP))
  rank_subset_overlap3[[i]] <- filter(aSNPs_haplotypes_v1_pap_fromfilt2_NEAND1KG, SNP %in% unique(rank_subsets[[i]]$SNP))
}


#proportions
rank_df <- data.frame(rank1 = c(0.0033,0.0054,0.0037), rank2 = c(0.0022,0.0056,0.0050), rank3 = c(0.0047,0.0041,0.010), rank4 = c(0.0031,0.0071,0.0071), rank5 = c(0.0055,0.0046,0.0073), ancestry = c("DEN", "NEANDPAP", "NEAND1KG"))

#change to long format
rank_long <- rank_df %>%
  pivot_longer(
    cols = starts_with("rank"),
    names_to = "rank",
    values_to = "proportion"
  )

#plot
pdf("aseQTL_archaic_prop.pdf", width=8, height=6)
ggplot(rank_long, aes(x = rank, y = proportion, color = ancestry, group = ancestry)) +
  geom_point(size = 3) +   
  geom_line() +             
  theme_bw() +        
  labs(
    x = "aseQTL rank",
    y = "Proportion of archaic aseQTLs"
  ) +  scale_color_manual(values = c("DEN" = "#2c6184","NEANDPAP" = "#e69b99","NEAND1KG" = "#89689d"))
  #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#test whether the observed variance in proportions across ranks is greater than what would be expected by chance
#perform a permutation test for each ancestry group
set.seed(123)

permutation_results <- rank_long %>%
  group_by(ancestry) %>%
  summarise(
    p_value = {
      observed_stat <- var(proportion)
      permuted_stats <- replicate(1000, var(sample(proportion)))
      mean(permuted_stats >= observed_stat)
    }
  )
print(permutation_results)



# ==========================
# Effect sizes of ASE SNPs
# ==========================

#subset to only include ASE SNPs
raw_cht2_gr_df2_filt2_mafs_ase <- raw_cht2_gr_df2_filt2_mafs[raw_cht2_gr_df2_filt2_mafs$Significant == "S",]
raw_cht2_gr_df2_filt2_mafs_ase$Effect <- abs(0.5 - raw_cht2_gr_df2_filt2_mafs_ase$RATIO)

#stats
stat.test <- raw_cht2_gr_df2_filt2_mafs_ase %>% t_test(Effect ~ All)
stat.test
stat.test <- stat.test %>% add_xy_position(x = "All")

#boxplots
pdf("ASE_effect_intro_vs_non.pdf")
bxp <- ggboxplot(raw_cht2_gr_df2_filt2_mafs_ase, x = "All", y = "Effect", fill = "All", 
                 palette = c("grey", "#E7B800", "#89689d"))
bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.shorten = 0.05)
dev.off()

