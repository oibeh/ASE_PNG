##############################################################
## Description: ASE + introgression (overlap and enrichment)
##############################################################


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
library(scales)
library(knitr)
library(karyoploteR)
library(regioneR)
library(zoo)
library(stringr)
library(purrr)
library(ggridges)
library(nullrangesData)
library(tidyr)
library(ExperimentHub)
library(AnnotationHub)
library(nullranges)
library(vcd)
library(corrplot)
library(forcats)
library(annotatr)


# ==================================
# Karyoplot for ASE + archaic SNPs
# ==================================

#get region dfs
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG[c("chr_st_end")]
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot <- data.frame(do.call('rbind', strsplit(as.character(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG$chr_st_end),'.',fixed=TRUE)))
colnames(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot) <- c("chrom","start","end")
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot %>% mutate_at(c('start', 'end'), as.integer)

aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN[c("chr_st_end")]
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot <- data.frame(do.call('rbind', strsplit(as.character(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN$chr_st_end),'.',fixed=TRUE)))
colnames(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot) <- c("chrom","start","end")
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot %>% mutate_at(c('start', 'end'), as.integer)

aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP[c("chr_st_end")]
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot <- data.frame(do.call('rbind', strsplit(as.character(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP$chr_st_end),'.',fixed=TRUE)))
colnames(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot) <- c("chrom","start","end")
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot %>% mutate_at(c('start', 'end'), as.integer)


#make granges objects
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot_gr <- GRanges(seqnames=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot$chrom,IRanges(start=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot$start,end=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot$end))
genome(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot_gr) <- 'hg38'
#why am I reducing here?? (for co-occurance analysis, didn't reduce)
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot_gr <- reduce(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot_gr)

aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot_gr <- GRanges(seqnames=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot$chrom,IRanges(start=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot$start,end=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot$end))
genome(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot_gr) <- 'hg38'
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot_gr <- reduce(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot_gr)

aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot_gr <- GRanges(seqnames=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot$chrom,IRanges(start=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot$start,end=aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot$end))
genome(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot_gr) <- 'hg38'
aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot_gr <- reduce(aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot_gr)


regions_NEAND1KG <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEAND1KG_annot_gr
regions_DEN <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_DEN_annot_gr
regions_NEANDPAP <- aSNPs_haplotypes_v1_pap_fromfilt_ASEcount_dstnct2_NEANDPAP_annot_gr

chroms <- str_sort(unique(all_introgression_calls_PNGx_annot_dstnct$seqnames), numeric=TRUE)
#c("Denisova" = "#2c6184", "Neand_1KG" = "#89689d", "Neand_PAP_uniq" = "#e69b99")

kp <- plotKaryotype(chromosomes="chr1")
kpPlotRegions(kp, data=regions_NEAND1KG, col="#89689d", border=NA, r0=0, r1=0.25)
kpPlotRegions(kp, data=regions_NEANDPAP, col="#e69b99", border=NA, r0=0.27, r1=0.52)
kpPlotRegions(kp, data=regions_DEN, col="#2c6184", border=NA, r0=0.54, r1=0.79)

#boxed regions, just introgression
pdf("newKP_region_allchroms_justintro.pdf", width=10)
kp <- plotKaryotype(genome="hg38",chromosomes=chroms)
kpPlotRegions(kp, data=regions_NEAND1KG, col="#89689d", border="#89689d", r0=0, r1=0.30)
kpPlotRegions(kp, data=regions_NEANDPAP, col="#e69b99", border="#e69b99", r0=0.31, r1=0.61)
kpPlotRegions(kp, data=regions_DEN, col="#2c6184", border="#2c6184", r0=0.62, r1=0.92)
dev.off()

#boxed regions with ASE density beneath
head(marks_pos) #these are the ASE SNPs
marks_pos_gr <- GRanges(seqnames=marks_pos$chrom,IRanges(start=marks_pos$start,width=1))
genome(marks_pos_gr) <- 'hg38'

pdf("newKP_region_allchroms_intro+ASE.pdf", width=10)
kp <- plotKaryotype(genome="hg38",chromosomes=chroms,plot.type=2)
kpPlotRegions(kp, data=regions_NEAND1KG, col="#89689d", border="#89689d", r0=0, r1=0.30)
kpPlotRegions(kp, data=regions_NEANDPAP, col="#e69b99", border="#e69b99", r0=0.31, r1=0.61)
kpPlotRegions(kp, data=regions_DEN, col="#2c6184", border="#2c6184", r0=0.62, r1=0.92)
kpPlotRegions(kp, data=marks_pos_gr,data.panel=2, col="darkgreen",r0=0, r1=0.3)
#kpPlotDensity(kp, data=marks_pos_gr,data.panel=2,r0=0, r1=0.8)
dev.off()


# ====================================
# Enrichment analysis with nullranges
# ====================================

#add sub_length and ID column to all_introgression_calls_PNGx_annot_gr_df
all_introgression_calls_PNGx_annot_gr_df$sub_length <- all_introgression_calls_PNGx_annot_gr_df$end - all_introgression_calls_PNGx_annot_gr_df$start
all_introgression_calls_PNGx_annot_gr_df$id <- seq(1:nrow(all_introgression_calls_PNGx_annot_gr_df))

#create Granges and keep lengths as metadata
library(BSgenome.Hsapiens.UCSC.hg38)

all_introgression_calls_PNGx_annot_granges <- GRanges(seqnames=all_introgression_calls_PNGx_annot_gr_df$seqnames,sub_length=all_introgression_calls_PNGx_annot_gr_df$sub_length, id=all_introgression_calls_PNGx_annot_gr_df$id, IRanges(start=all_introgression_calls_PNGx_annot_gr_df$start,end=all_introgression_calls_PNGx_annot_gr_df$end), seqlengths=seqlengths(Hsapiens))
genome(all_introgression_calls_PNGx_annot_granges) <- 'hg38'
length(all_introgression_calls_PNGx_annot_granges)
table(seqnames(all_introgression_calls_PNGx_annot_granges))

#use pre-built hg38 segmentation by gene density for hg38 build
eh <- ExperimentHub()
seg_cbs <- eh[["EH7307"]] #prefer CBS for hg38
seg_hmm <- eh[["EH7308"]] #HMM segmentation
seg <- seg_cbs


set.seed(5) 
seg <- eh[["EH7307"]] #pre-built genome segmentation for hg38
#input `ranges` require seqlengths, if missing see `GenomeInfoDb::Seqinfo`
seqlengths(all_introgression_calls_PNGx_annot_granges)
#next, we remove non-standard chromosomes...
all_introgression_calls_PNGx_annot_granges <- keepStandardChromosomes(all_introgression_calls_PNGx_annot_granges, pruning.mode="coarse")
#...and mitochondrial genome + X, Y, as these are too short/not needed
seqlevels(all_introgression_calls_PNGx_annot_granges, pruning.mode="coarse") <- setdiff(seqlevels(all_introgression_calls_PNGx_annot_granges), c("chrM", "MT"))

R <- 1000
blockLength <- 5.5e6 #use 5.5 for now....optimal
boots_mine <- bootRanges(all_introgression_calls_PNGx_annot_granges, proportionLength=TRUE, blockLength, R = R, seg = seg, exclude=exclude)
boots_mine

#(Note that in this case blockLength is the maximal block length possible, but the actual block lengths per segment 
#will be smaller---proportional to the fraction of basepairs covered by that state in the genome segmentation

#assessing quality of bootstrap samples
#We can examine properties of permuted y over iterations, and compare to the original y. 
#To do so, we first add the original features as iter=0. Then compute summaries:
combined <- all_introgression_calls_PNGx_annot_granges %>% 
  mutate(iter=0) %>%
  bind_ranges(boots_mine) %>% 
  plyranges::select(iter)
stats <- combined %>% 
  group_by(iter) %>%
  summarize(n = n()) %>%
  as_tibble()
head(stats) #the number of ranges in the granges objects (iter 0 is the original observed data that we just added)

# We can also look at distributions of various aspects, 
#e.g. here the inter-feature distance of features, across a few of the bootstraps and the original feature set y.
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
interdist <- function(dat) {
    x <- dat[-1,]
    y <- dat[-nrow(dat),]
    ifelse(x$seqnames == y$seqnames,
           x$start + floor((x$width - 1)/2) -
           y$start - floor((y$width - 1)/2), NA)
}

#just looking at first 15 iterations...
pdf("nullranges_interdist_15.pdf")
combined %>% plyranges::filter(iter %in% 0:15) %>%
  mutate(iter = droplevels(iter)) %>%
  plyranges::select(iter) %>%
  as_tibble() %>% 
  nest(data = !iter) %>%
  mutate(interdist = map(data, interdist)) %>% 
  dplyr::select(iter, interdist) %>% 
  unnest(interdist) %>% 
  mutate(type = ifelse(iter == 0, "original", "boot"),
         interdist = pmax(interdist, 0)) %>%
  filter(!is.na(interdist)) %>%
  ggplot(aes(log10(interdist + 1), iter, fill=type)) +
  geom_density_ridges(alpha = 0.75) +
  geom_text(data = head(stats, 16),
            aes(x=1.5, y=iter, label=paste0("n=",n), fill=NULL),
            vjust=1.5) + theme_bw() + xlab(expression('Log'[10]*'(interdist + 1)')) + ylab("Iteration")


#Now we can combine bootstrapping with plyranges to perform statistical enrichment analysis
#The general idea will be to combine the long vector of bootstrapped ranges, indexed by iter, 
#with another set of ranges to compute enrichment
boot_stats <- raw_cht_padjpass_gr_df_pos_granges %>% join_overlap_inner(boots_mine) %>%
  group_by(x_id, iter) %>%
  summarize(n_overlaps = n()) %>%
  as_tibble() %>%
  complete(x_id, iter, fill=list(n_overlaps = 0)) %>%
  group_by(iter) %>%
  summarize(sumOverlaps = sum(n_overlaps))

#plot a histogram
pdf("nullranges_bootstats_distr.pdf")
ggplot(boot_stats, aes(sumOverlaps)) +
  geom_histogram(binwidth=5)+
  geom_vline(xintercept = sum(raw_cht_padjpass_gr_df_pos_granges$n_overlaps), linetype = "dashed", col="blue") + theme_bw() + xlab("No. of introgressed ASE SNPs") + ylab("Count")
dev.off()


# ==========================================================
# Contrasting ASE at introgressed and non-introgressed SNPs
# ==========================================================

#snps with mafs
p1.snp2 <- as_granges(snps_selection, seqnames=TEST.SNP.CHROM, start=TEST.SNP.POS, end=TEST.SNP.POS)

#find intersect for each ancestry
overlaps_selection_NEAND1KG <- join_overlap_inner(p1.snp2, regions_NEAND1KG)
overlaps_selection_NEAND1KG <- as.data.frame(overlaps_selection_NEAND1KG)
dim(overlaps_selection_NEAND1KG)
overlaps_selection_NEAND1KG$SNP <- paste(overlaps_selection_NEAND1KG$seqnames, overlaps_selection_NEAND1KG$start, sep=":")
overlaps_selection_NEAND1KG$NEAND <- rep("NEAND1KG", length(overlaps_selection_NEAND1KG$SNP))

overlaps_selection_NEANDPAP <- join_overlap_inner(p1.snp2, regions_NEANDPAP)
overlaps_selection_NEANDPAP <- as.data.frame(overlaps_selection_NEANDPAP)
dim(overlaps_selection_NEANDPAP)
overlaps_selection_NEANDPAP$SNP <- paste(overlaps_selection_NEANDPAP$seqnames, overlaps_selection_NEANDPAP$start, sep=":")
overlaps_selection_NEANDPAP$NEANDPAP <- rep("NEANDPAP", length(overlaps_selection_NEANDPAP$SNP))

overlaps_selection_DEN <- join_overlap_inner(p1.snp2, regions_DEN)
overlaps_selection_DEN <- as.data.frame(overlaps_selection_DEN)
dim(overlaps_selection_DEN)
overlaps_selection_DEN$SNP <- paste(overlaps_selection_DEN$seqnames, overlaps_selection_DEN$start, sep=":")
overlaps_selection_DEN$DEN <- rep("DEN", length(overlaps_selection_DEN$SNP))

#merge with full ASE stats
raw_cht2_gr_df2_filt2_mafs$NEAND <- overlaps_selection_NEAND1KG$NEAND[match(raw_cht2_gr_df2_filt2_mafs$SNP, overlaps_selection_NEAND1KG$SNP)]
raw_cht2_gr_df2_filt2_mafs$NEANDPAP <- overlaps_selection_NEANDPAP$NEANDPAP[match(raw_cht2_gr_df2_filt2_mafs$SNP, overlaps_selection_NEANDPAP$SNP)]
raw_cht2_gr_df2_filt2_mafs$DEN <- overlaps_selection_DEN$DEN[match(raw_cht2_gr_df2_filt2_mafs$SNP, overlaps_selection_DEN$SNP)]
raw_cht2_gr_df2_filt2_mafs$NEAND[is.na(raw_cht2_gr_df2_filt2_mafs$NEAND)] <- "No"
raw_cht2_gr_df2_filt2_mafs$NEANDPAP[is.na(raw_cht2_gr_df2_filt2_mafs$NEANDPAP)] <- "No"
raw_cht2_gr_df2_filt2_mafs$DEN[is.na(raw_cht2_gr_df2_filt2_mafs$DEN)] <- "No"

#make a column for both DEN sources together
raw_cht2_gr_df2_filt2_mafs$DEN_both <- ifelse(raw_cht2_gr_df2_filt2_mafs$NEANDPAP == "NEANDPAP" | raw_cht2_gr_df2_filt2_mafs$DEN == "DEN", "D", "No")

#run for Den then Neand
propDiffBin <- function(dt, minFreq, maxFreq) {
	dtSubset <- dt[dt$maf >= minFreq & dt$maf < maxFreq,]
	n_sig_case <- nrow(dtSubset[dtSubset$DEN_both == "D" & dtSubset$Significant == "S",])
	n_ns_case <- nrow(dtSubset[dtSubset$DEN_both == "D" & dtSubset$Significant == "NS",])
	n_sig_control <- nrow(dtSubset[dtSubset$DEN_both == "No" & dtSubset$Significant == "S",])
	n_ns_control <- nrow(dtSubset[dtSubset$DEN_both == "No" & dtSubset$Significant == "NS",])
	S <- 100000
	N <- rbeta(S, n_sig_case + 1, n_ns_case + 1)
	H <- rbeta(S, n_sig_control + 1, n_ns_control + 1)
	interval <- quantile(N - H, c(0.025, 0.5, 0.975))
	return(data.table(min = minFreq, max = maxFreq, 
	n_introgressed = nrow(dtSubset[dtSubset$DEN_both == "D",]), 
	n_nonintrogressed = nrow(dtSubset[dtSubset$DEN_both == "No",]), 
	CI0.025 = interval[1], CI0.5 = interval[2], CI0.975 = interval[3]))
}

diffByDAF <- do.call(rbind, lapply(seq(0, 0.5, 0.1), function(x) propDiffBin(raw_cht2_gr_df2_filt2_mafs, x, x + 0.1)))
limits <- aes(ymin = CI0.025, ymax = CI0.975, x = (min + max) / 2)
#plot
pdf("MAF_prop_1.pdf")
ggplot(data = diffByDAF[n_introgressed >= 20], aes(x = (min + max) / 2, y = CI0.5)) + 
	geom_point() + 
	geom_errorbar(limits) +
	theme_bw() +
	ylab(expression(paste(italic("p")["Den"] - italic("p")["M"], " (95 % CI)"))) +
	xlab("Minor allele frequency (PAP)") +
	geom_hline(yintercept = 0, lty = "dashed", color = "red")
dev.off()

