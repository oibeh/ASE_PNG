#########################################
## Description: PNG introgression calls
#########################################

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


# =============================================
# Load Danat's introgression files and filter
# =============================================

aSNPs_haplotypes_v1 <- read.table("/Users/nibeh/Documents/PNG_ASE/From_Danat/aSNPs.haplotypes.v1.tsv",header=T,check.names=FALSE)
all_introgression_calls <- read.table("/Users/nibeh/Documents/PNG_ASE/From_Danat/Isabela/PAP_10_EUR_7.txt",header=T,check.names=FALSE) #Aug 19th, all samples, not filtered
head(aSNPs_haplotypes_v1)

#grab PAP rows
aSNPs_haplotypes_v1_pap <- aSNPs_haplotypes_v1[aSNPs_haplotypes_v1$SupPop == "PAP",]
colnames(aSNPs_haplotypes_v1_pap)[20] <- "SupPop2" #rename so it doesn't interfere with the SupPop...not sure if needed but will keep doing this

#now filter the haplotypes from Isabela (initially she filtered out Nean, but I'm keeping all)
head(all_introgression_calls)
dim(all_introgression_calls) #25572 (aSNPs in all populations)
all_introgression_calls <- all_introgression_calls[all_introgression_calls$SupPop == "PAP",]
dim(all_introgression_calls) #19033
#add Melbourne_ID to missing rows
#df$spring[grepl("W02$",df$Time)] <- 1
all_introgression_calls$Melbourne_ID[grepl("PNG8",all_introgression_calls$haplo)] <- "PNG8"
all_introgression_calls$Melbourne_ID[grepl("PNG10",all_introgression_calls$haplo)] <- "PNG10"
all_introgression_calls$Melbourne_ID[grepl("PNG16",all_introgression_calls$haplo)] <- "PNG16"
all_introgression_calls$Melbourne_ID[grepl("PNG25",all_introgression_calls$haplo)] <- "PNG25"
all_introgression_calls$Melbourne_ID[grepl("PNG26",all_introgression_calls$haplo)] <- "PNG26"
#remove png21
all_introgression_calls <- all_introgression_calls[!grepl('PNG21', all_introgression_calls$Melbourne_ID),]


#Next, get all introgressed regions across samples, reduce, and get gene names
head(all_introgression_calls)

#use the set of all samples
all_introgression_calls_PNGx_annot <- data.frame(do.call('rbind', strsplit(as.character(all_introgression_calls$chr_st_end),'.',fixed=TRUE)))
#reformat the dfs
colnames(all_introgression_calls_PNGx_annot) <- c("chrom","start","end")
#fix class
str(all_introgression_calls_PNGx_annot)
all_introgression_calls_PNGx_annot <- all_introgression_calls_PNGx_annot %>% mutate_at(c('start', 'end'), as.integer)
#check
str(all_introgression_calls_PNGx_annot)

#make a GR object with the overlapping regions (instead of doing it one by one)
#merging ranges across individuals
dim(all_introgression_calls_PNGx_annot)

#create granges, reduce, and then annotate to get gene names
all_introgression_calls_PNGx_annot_gr <- GRanges(seqnames=all_introgression_calls_PNGx_annot$chrom,IRanges(start=all_introgression_calls_PNGx_annot$start,end=all_introgression_calls_PNGx_annot$end))
genome(all_introgression_calls_PNGx_annot_gr) <- 'hg38'
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#reduce
all_introgression_calls_PNGx_annot_gr <- reduce(all_introgression_calls_PNGx_annot_gr)
peakAnno <- annotatePeak(all_introgression_calls_PNGx_annot_gr, tssRegion=c(-1000, 1000), TxDb=txdb, level="gene", annoDb="org.Hs.eg.db",overlap="all", verbose=FALSE)
length(as.data.frame(peakAnno)$SYMBOL)
all_introgression_calls_PNGx_annot_gr_df <- as.data.frame(all_introgression_calls_PNGx_annot_gr)
all_introgression_calls_PNGx_annot_gr_df <- all_introgression_calls_PNGx_annot_gr_df[,1:3]
all_introgression_calls_PNGx_annot_gr_df$SYMBOL <- as.data.frame(peakAnno)$SYMBOL


#filter for 10 aSNPs per haplotype
all_introgression_calls_filt <- all_introgression_calls[all_introgression_calls$aSNPs_num >= 10,]
dim(all_introgression_calls_filt)

#get aSNP rows in aSNPs_haplotypes_v1_pap that are found in all_introgression_calls_filt using chr_st_end column
aSNPs_haplotypes_v1_pap_fromfilt <- filter(aSNPs_haplotypes_v1_pap, chr_st_end %in% unique(all_introgression_calls_filt$chr_st_end))



