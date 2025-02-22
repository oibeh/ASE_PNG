###################################################
## Description: PNG ASE processing and annotation
###################################################


# =================
# Load libraries
# =================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(colourvalues)
library(RColorBrewer)
library(ggExtra)
library(rstatix)
library(patchwork)
library(tidyverse)
library(valr)
library(ChIPseeker)
library(GenomicRanges)
library(plyranges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(UpSetR)
library(devtools)
library(ggpubr)
library(cowplot)
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


# ====================================
# Summarize and annotate WASP output
# ====================================

#read in raw CHT results
load("/Users/nibeh/Documents/PNG_ASE/WASP/CHT_full/New_exonregions/005_v2/WASP005v2.Rdata")
raw_cht <- read.table("/Users/nibeh/Documents/PNG_ASE/WASP/CHT_full/New_exonregions/cht_results.txt", header=TRUE)

head(raw_cht)
#add SNP column for merging downstream
raw_cht$SNP <- paste(raw_cht$TEST.SNP.CHROM, raw_cht$TEST.SNP.POS, sep=":")
#remove rows with no rsID
raw_cht <- raw_cht[raw_cht$TEST.SNP.ID != "b'.'", ]

#FDR adjustment
raw_cht <- adjust_pvalue(raw_cht, p.col = "P.VALUE", output.col = "P.ADJ", method = "BH")
#how many pass?
raw_cht_pvalpass <- raw_cht[raw_cht$P.VALUE <= 0.05,]
raw_cht_padjpass <- raw_cht[raw_cht$P.ADJ <= 0.05,] #using 0.05 now
dim(raw_cht_padjpass)
#unique SNPs
length(unique(raw_cht_padjpass$TEST.SNP.ID))

#next, map SNPs to genes
raw_cht_padjpass_gr <- GRanges(seqnames=raw_cht_padjpass$TEST.SNP.CHROM,IRanges(start=raw_cht_padjpass$TEST.SNP.POS,width=1))
genome(raw_cht_padjpass_gr) <- 'hg38'
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(raw_cht_padjpass_gr, tssRegion=c(-1000, 1000), TxDb=txdb, level="gene", annoDb="org.Hs.eg.db",overlap="all", verbose=FALSE)
length(as.data.frame(peakAnno)$SYMBOL)
raw_cht_padjpass_gr_df <- as.data.frame(raw_cht_padjpass_gr)
raw_cht_padjpass_gr_df <- raw_cht_padjpass_gr_df[,1:2]
raw_cht_padjpass_gr_df$Symbol <- as.data.frame(peakAnno)$SYMBOL
ASE_genes <- unique(raw_cht_padjpass_gr_df$Symbol)
length(ASE_genes)

#create SNP column for merging afterwards
raw_cht_padjpass_gr_df$SNP <- paste(raw_cht_padjpass_gr_df$seqnames, raw_cht_padjpass_gr_df$start, sep=":")
colnames(raw_cht_padjpass_gr_df) <- c("Chr","Start","Symbol","SNP")
#add genes to main df with sig ASE SNPs
raw_cht_padjpass$SYMBOL <- raw_cht_padjpass_gr_df$Symbol
raw_cht_padjpass_gr_df2 <- raw_cht_padjpass #duplicate


# =================
# ASE QC plots (I)
# =================

#Plotting reference allele mapping ratio
#the goal is to check if we have overcome the intrinsic reference allele bias observed in ASE testing
raw_cht2 <- raw_cht
raw_cht2$RATIO <- raw_cht2$REF.AS.READ.COUNT/raw_cht2$TOTAL.AS.READ.COUNT
p <- ggplot(raw_cht2, aes(x=RATIO)) + geom_density(fill="grey") + geom_vline(aes(xintercept=mean(RATIO)),
            color="#21918c", linetype="dashed", linewidth=1) + theme_bw(base_size = 18) + ylab("Density") + xlab("Reference allele mapping ratio") + theme_bw()

pdf("Ref_allele_ratio.pdf")
p
dev.off()


#Plotting mapping ratio vs. total coverage
#first, add gene names to full table (raw_cht2)
head(raw_cht2)

raw_cht2_gr <- GRanges(seqnames=raw_cht2$TEST.SNP.CHROM,IRanges(start=raw_cht2$TEST.SNP.POS,width=1))
genome(raw_cht2_gr) <- 'hg38'
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(raw_cht2_gr, tssRegion=c(-1000, 1000), TxDb=txdb, level="gene", annoDb="org.Hs.eg.db",overlap="all", verbose=FALSE)
length(as.data.frame(peakAnno)$SYMBOL)
raw_cht2_gr_df <- as.data.frame(raw_cht2_gr)
raw_cht2_gr_df <- raw_cht2_gr_df[,1:2]
raw_cht2_gr_df$Symbol <- as.data.frame(peakAnno)$SYMBOL
raw_cht2_gr_df$SNP <- paste(raw_cht2_gr_df$seqnames, raw_cht2_gr_df$start, sep=":")

raw_cht2_gr_df2 <- raw_cht2
raw_cht2_gr_df2$SYMBOL <- raw_cht2_gr_df$Symbol

#now add significance for each SNP
raw_cht2_gr_df2 <- raw_cht2_gr_df2 %>% mutate(Significant = if_else(P.ADJ <= 0.05, "S", "NS")) #changed to 0.05

p2 <- ggplot(raw_cht2_gr_df2,aes(x=RATIO,y=TOTAL.AS.READ.COUNT, col=Significant))+geom_point() + 
theme_bw(base_size = 18) + ylab("Coverage") + xlab("Reference allele mapping ratio") + scale_y_log10()

#change ordering of point being plotted
#with genes labeled
p2 <- ggplot(raw_cht2_gr_df2 %>% arrange(Significant), label=SYMBOL) + geom_point(aes(x=RATIO,y=TOTAL.AS.READ.COUNT, col=Significant)) + 
theme_bw(base_size = 18) + ylab("Coverage") + xlab("Reference allele mapping ratio") + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + scale_colour_manual("ASE Significance", values = c("grey", "#21918c")) +
geom_text_repel(data=subset(raw_cht2_gr_df2, TOTAL.AS.READ.COUNT > 80000), aes(RATIO, TOTAL.AS.READ.COUNT, label = SYMBOL), show.legend = FALSE, nudge_y=0.9, verbose=FALSE, max.overlaps = 150)

pdf("ASE_significance_withgenes.pdf", width=14)
p2
dev.off()

#without genes labeled
p2 <- ggplot(raw_cht2_gr_df2 %>% arrange(Significant), label=Symbol) + geom_point(aes(x=RATIO,y=TOTAL.AS.READ.COUNT, col=Significant)) + 
theme_bw(base_size = 18) + ylab("Coverage") + xlab("Reference allele mapping ratio") + scale_y_log10(labels = function(x) format(x, scientific = FALSE), limits = c(80,600000)) + scale_colour_manual("ASE Significance", values = c("grey", "#21918c"))

pdf("ASE_significance_withoutgenes.pdf", width=14)
p2
dev.off()


#Cumulative distribution of ASE SNPs
head(table(raw_cht_padjpass$SYMBOL))
gene_tab <- data.frame(cbind(table(raw_cht_padjpass$SYMBOL)))

pdf("ASE_gene_histogram.pdf")
ggplot(gene_tab, aes(x = gene_tab)) + geom_histogram(colour = 1, fill = "white", bins = 200) + theme_bw() +
ylab("Number of genes") + xlab("Number of ASE SNPs") + scale_y_log10()
dev.off()


#Mean TPM against SNP coverage
qc2 <- ggplot(raw_cht2_gr_df2_filt, aes(x=Mean_TPM, y=TOTAL.AS.READ.COUNT)) + 
    geom_point(alpha = .1) + geom_smooth() + theme_bw() + ylab("Coverage") + xlab(bquote(Log[2](TPM)))
suppressWarnings(print(qc2))

#log10 y
qc3 <- ggplot(raw_cht2_gr_df2_filt, aes(x=Mean_TPM, y=TOTAL.AS.READ.COUNT)) + 
    geom_point(alpha = .1) + geom_smooth() + theme_bw() + ylab("Coverage") + xlab(bquote(Log[2](TPM))) + scale_y_continuous(trans='log10') 
suppressWarnings(print(qc3))

#plot significance (cut-off is marked) against effect size 
qc4 <- ggplot(raw_cht2_gr_df2_edit[raw_cht2_gr_df2_edit$NegLog10P < 200,], aes(x=Effect, y=NegLog10P)) + 
    geom_point(alpha = .1) + theme_bw() + xlab("Effect (deviation from 0.5 allele balance)") + ylab(bquote(-Log[10](p-value)))+
    geom_hline(yintercept = -log10(max(raw_cht_padjpass_gr_df2$P.VALUE)), linetype = "dashed", color = "darkblue", size = 0.7)
suppressWarnings(print(qc4))

#plot significance (cut-off is marked) against expression
qc5 <- ggplot(raw_cht2_gr_df2_edit[raw_cht2_gr_df2_edit$NegLog10P < 200,], aes(x=Mean_TPM, y=NegLog10P)) + 
    geom_point(alpha = .1) + theme_bw() + xlab(bquote(Log[2](TPM))) + ylab(bquote(-Log[10](p-value)))
suppressWarnings(print(qc5))

#proportion of ASE SNPs in ASE genes across TPM bins
level_ord <- levels(bins)
p_prop <- propASE_wTPM %>% 
  ggplot2::ggplot(., aes(y = Proportion_S, x = factor(bins, levels = level_ord))) +
  geom_point(colour= "#0000FF",alpha = 0.2, position = position_jitter(width = 0.1, seed = 3922)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot") + 
  stat_summary(aes(y = Proportion_S, group = 1),fun.data = "mean_cl_boot", geom="line", size=1.1, show.legend=FALSE) +
  theme_bw() + ylab("Proportion of ASE SNPs in ASE genes") + xlab(bquote(Log[2](TPM)))
suppressWarnings(print(p_prop))



# =================
# ASE QC plots (II)
# =================

#Number of ASE SNPs vs ASE gene length

#get exon annotations
#connect to the Ensembl database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#retrieve exon annotations for the ASE genes
exon_annotations <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end",
                 "external_gene_name", "strand"),
  filters = "external_gene_name",
  values = ase_genes,
  mart = ensembl
)
#convert to GRanges object
exon_granges <- GRanges(
  seqnames = paste0("chr", exon_annotations$chromosome_name),
  ranges = IRanges(start = exon_annotations$exon_chrom_start, end = exon_annotations$exon_chrom_end),
  strand = exon_annotations$strand,
  gene_name = exon_annotations$external_gene_name
)

#split exons by gene name
exons_split <- split(exons_of_interest, exons_of_interest$gene_name)
#get union of exons
reduced_exons <- lapply(exons_split, reduce)
#calculate the total length for each gene
gene_lengths <- sapply(reduced_exons, function(exon) sum(width(exon)))
gene_length_df <- data.frame(Gene = names(gene_lengths), Length = gene_lengths)
print(gene_length_df)

#get ASE SNP counts for each gene
raw_cht2_gr_df2_edit2_sig_un <- raw_cht2_gr_df2_edit2_sig[!duplicated(raw_cht2_gr_df2_edit2_sig$SNP), ]
SNP_counts <- as.data.frame(table(raw_cht2_gr_df2_edit2_sig_un$SYMBOL))
colnames(SNP_counts) <- c("hgnc_symbol","ASEsnp_count")
#merge
colnames(gene_length_df) <- c("hgnc_symbol","Length")
SNP_counts$Length <- gene_length_df$Length[match(SNP_counts$hgnc_symbol,gene_length_df$hgnc_symbol)]
SNP_counts <- SNP_counts[order(SNP_counts$ASEsnp_count, decreasing = TRUE), ]
SNP_counts$label <- c(rep("top50",50), rep("rest", 2034))

pdf("ASE_genelength.pdf")
ggplot(SNP_counts, aes(x=Length, y=ASEsnp_count, color = label)) + 
    geom_point(alpha = .3) + scale_x_continuous(trans='log10') + theme_bw() + labs(color = "Category") +
    scale_color_manual(values = c("top50" = "#F8766D", "rest" = "#333333")) +
    geom_label_repel(data = SNP_counts[SNP_counts$label == "top50", ], 
                   aes(label = hgnc_symbol), 
                   box.padding = 0.5, 
                   point.padding = 0.5,  
                   segment.color = '#F8766D',  
                   show.legend = FALSE, max.overlaps =40) + xlab("Gene length") + ylab("Number of ASE SNPs")
dev.off()

#now add proportion counts to ASEsnp_count
gene_counts <- raw_cht2_gr_df2_filt %>%
  count(SYMBOL, name = "row_count")

SNP_counts2 <- SNP_counts2 %>%
  left_join(gene_counts, by = c("hgnc_symbol" = "SYMBOL"))

#replace NA values in the new column with 0 for genes with no match
SNP_counts2 <- SNP_counts2 %>%
  mutate(row_count = ifelse(is.na(row_count), 0, row_count))

SNP_counts2$prop <- (SNP_counts2$ASEsnp_count)/(SNP_counts2$row_count)

pdf("ASE_genelength_prop.pdf")
ggplot(SNP_counts2, aes(x=size, y=prop, color = label)) + 
    geom_point(alpha = .3) + scale_x_continuous(trans='log10') + theme_bw() + labs(color = "Category") +
    scale_color_manual(values = c("top50" = "#F8766D", "rest" = "#333333")) +
    geom_label_repel(data = SNP_counts2[SNP_counts2$label == "top50", ], 
                   aes(label = hgnc_symbol), 
                   box.padding = 0.5,  
                   point.padding = 0.5,  
                   segment.color = '#F8766D',  
                   show.legend = FALSE, max.overlaps =40) + xlab("Gene length") + ylab("Proportion of signficant ASE SNPs")
dev.off()



# ===================
# ASE QC plots (III)
# ===================

#ASE concordance analysis (exon level)

p_values <- 10^seq(log10(0.01), log10(5e-9), by = -1)
#initialize an empty vector to store the discordance rates
discordance_rates <- numeric(length(p_values))

for (i in seq_along(p_values)) {

	p_value <- p_values[i]

	head(raw_cht2_gr_df2_edit2)
	#take all rows with pval
	raw_cht2_gr_df2_edit2_sig <- raw_cht2_gr_df2_edit2[raw_cht2_gr_df2_edit2$P.VALUE < p_value,]

	#remove genes that are present only once in the df
	raw_cht2_gr_df2_edit2_sig <- raw_cht2_gr_df2_edit2_sig %>% group_by(SYMBOL) %>% dplyr::filter(n() > 1) %>% ungroup()
	raw_cht2_gr_df2_edit2_sig <- as.data.frame(raw_cht2_gr_df2_edit2_sig)
	genes_unique <- unique(raw_cht2_gr_df2_edit2_sig$SYMBOL)

	gg <- annotLookup[annotLookup$hgnc_symbol %in% genes_unique,]
	gg <- gg[!duplicated(gg$hgnc_symbol), ]
	entrez_ids_ase <- gg$entrezgene_id

	#subset exon list using entrez_ids_ase
	x_filt <- x[x$GeneID %in% entrez_ids_ase,]

	#find which exons the ASE SNPs fall into
	raw_cht2_gr_df2_edit2_sig_unique <- raw_cht2_gr_df2_edit2_sig[!duplicated(raw_cht2_gr_df2_edit2_sig$SNP), ]
	dim(raw_cht2_gr_df2_edit2_sig_unique)

	#actually, don't need granges. get coordinate and use valr instead
	x_filt_coord <- x_filt[c("Chr","Start","End")]
	colnames(x_filt_coord) <- c("chrom","start","end")

	snps_ase_pc_coord <- raw_cht2_gr_df2_edit2_sig_unique[c("TEST.SNP.CHROM","TEST.SNP.POS","TEST.SNP.POS")]
	colnames(snps_ase_pc_coord) <- c("chrom","start","end")

	exons_intersect <- valr::bed_intersect(x_filt_coord, snps_ase_pc_coord, suffix = c("_1", "_2"))
	exons_intersect <- as.data.frame(exons_intersect)

	#add more info
	exons_intersect$chr_st_end1 <- paste(exons_intersect$chrom, exons_intersect$start_1, exons_intersect$end_1, sep=".")
	exons_intersect$chr_st_end2 <- paste(exons_intersect$chrom, exons_intersect$start_2, exons_intersect$end_2, sep=".")
	exons_intersect$SNP <- paste(exons_intersect$chrom, exons_intersect$start_2, sep=":")

	#now add RATIO to exons_intersect
	exons_intersect$RATIO <- raw_cht2_gr_df2_edit2_sig_unique$RATIO[match(exons_intersect$SNP, raw_cht2_gr_df2_edit2_sig_unique$SNP)]

	#ok, now run the concordance test using exons_intersect
	#remove exons that occur only once
	exon_counts <- table(exons_intersect$chr_st_end1)
	exons_intersect <- exons_intersect[exons_intersect$chr_st_end1 %in% names(exon_counts[exon_counts > 1]), ]

	#get two columns
	exons_intersect_sub <- exons_intersect[c("RATIO","chr_st_end1")]

	data_exons <- exons_intersect_sub %>%
  	group_by(chr_st_end1) %>%
  	summarise(RATIO = t(combn(RATIO, 2)))
	out_exons <- cbind(data_exons[1], data.frame(data_exons$RATIO))
	out_exons <- as.data.frame(out_exons)

	out_exons2 <- out_exons %>%
   		mutate(label = case_when(
      	(X1 < 0.5 & X2 > 0.5) | (X1 > 0.5 & X2 < 0.5) ~ "grey",
      	(X1 <= 0.5 & X2 <= 0.5) | (X1 >= 0.5 & X2 >= 0.5) ~ "darkgreen"
	))

	rate_discord <- sum(out_exons2$label == "grey") / nrow(out_exons2)

	discordance_rates[i] <- rate_discord
}




