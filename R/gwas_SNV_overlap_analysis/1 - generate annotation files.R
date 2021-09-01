#
#         _       _    _   _         _                    _             _        
#        / /\    / /\ /\_\/\_\ _    / /\                /\ \           /\_\      
#       / / /   / / // / / / //\_\ / /  \              /  \ \         / / /  _   
#      / /_/   / / //\ \/ \ \/ / // / /\ \            / /\ \ \       / / /  /\_\ 
#     / /\ \__/ / //  \____\__/ // / /\ \ \          / / /\ \ \     / / /__/ / / 
#    / /\ \___\/ // /\/________// / /  \ \ \        / / /  \ \_\   / /\_____/ /  
#   / / /\/___/ // / /\/_// / // / /___/ /\ \      / / /    \/_/  / /\_______/   
#  / / /   / / // / /    / / // / /_____/ /\ \    / / /          / / /\ \ \      
# / / /   / / // / /    / / // /_________/\ \ \  / / /________  / / /  \ \ \     
#/ / /   / / / \/_/    / / // / /_       __\ \_\/ / /_________\/ / /    \ \ \    
#\/_/    \/_/          \/_/ \_\___\     /____/_/\/____________/\/_/      \_\_\   
# Script by Harry MacKay 
# hamackay@bcm.edu
# This script is the final revision of the single-SNP matching script
# It's functionally the same as the one in the first folder, but I've re-written it to be cleaner

library(tidyverse)
library(data.table)
library(ggplot2)
library(stringr)
library(corpora)

# The following establishes the HM Signature Graph Style (tm)
theme_hm <- function (base_size = 14) 
{
  half_line <- base_size/2
  theme_grey(base_size = base_size) %+replace% 
    theme(axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.text = element_text(size = rel(0.8), colour = "black"),
          axis.title.x = element_text(size=base_size-3, margin = margin(t = 1 * half_line, b = 0.8 * half_line/2)),
          axis.title.y = element_text(size=base_size-3, angle = 90, margin = margin(r = 1.5 * half_line, l = 0.8 * half_line/2)),
          axis.ticks = element_line(colour = "black"), 
          axis.ticks.length = unit(half_line/2, "pt"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_blank(), 
          #panel.grid.major = element_blank(),
          #panel.grid.major = element_line(colour = "grey92", size=0.5), 
          #panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
          strip.background = element_rect(fill = "grey85", colour = "grey20"), 
          legend.key = element_rect(fill = "white", colour = NA), 
          axis.text.y=element_text(size=base_size-2),
          axis.text.x=element_text(size=base_size-2),
          complete = TRUE)
}

# This is an all-purpose saving function that can be called upon to save the present graph in the correct format and folder
func_ggsave <- function(plot_name, gen_name, w_dim, h_dim){
  png_name <- paste("graph_png/", gen_name, ".png", sep="")
  pdf_name <- paste("graph_pdf/", gen_name, ".pdf", sep="")
  ggsave(plot_name, filename = png_name, width = w_dim, height = h_dim, dpi = 300, units = "in")
  ggsave(plot_name, filename = pdf_name, width = w_dim, height = h_dim, dpi = 300, units = "in")
}

func_legendsave <- function(plot_name, gen_name, w_dim, h_dim){
  png_name <- paste("legend_png/", gen_name, ".png", sep="")
  pdf_name <- paste("legend_pdf/", gen_name, ".pdf", sep="")
  ggsave(plot_name, filename = png_name, width = w_dim, height = h_dim, dpi = 300, units = "in")
  ggsave(plot_name, filename = pdf_name, width = w_dim, height = h_dim, dpi = 300, units = "in")
}

options(scipen = 999) # This is essential.  It prevents R from outputting genomic coordinates in useless scientific notation
# To make this easier to run on multiple platforms, I've made the directories variable
path_prefix <- "S:/" # Workstation
#path_prefix <- "/Volumes/bcm-pedi-nutrition-waterlandlab$/" # Macbook
#path_prefix <- "/mnt/share/" # Waterland2

setwd(paste(path_prefix, "DATA/Harry/Coarfa SNP project/NHGRI Enrichment rev4", sep=""))

# Generate a list of chromosomes
chrom_list <- seq(from=1, to=22, by=1)
chrom_ids <- paste0("chr", chrom_list)

# A list of tissue types
tissue_list <- c("blood", "brain" , "lung", "nerve", "skin", "thyroid")

# How big a window do you want around the corsiv?
# We're using 1mb windows for the final analysis since those match the ones Coarfa set up
corsiv_window <- 1000000

# Load the frequency file from the 1000 genomes project (but keeping just the variants that were tested)
freq_df <- fread(paste0(path_prefix, "DATA/Harry/Coarfa SNP project/NHGRI Enrichment rev2/output_csv/combined allele frequency - gtex selection.csv"))
freq_df$A1_freq <- 1-freq_df$maf
freq_df$A2_freq <- freq_df$maf

# It's useful to have the chromosome lengths since sometimes the act of expanding a window around the corsivs extends beyond the length of their chromosome
hg38_chrom_lengths <- read_csv((paste0(path_prefix, "DATA/Harry/Coarfa SNP project/NHGRI Enrichment rev2/hg38_chrom_lengths.csv")))
hg38_chrom_lengths <- subset(hg38_chrom_lengths, chr!="chrX" & chr!="chrY")

# Load the latest NHGRI trait SNPs
gwas_cat <- readr::read_delim((paste0(path_prefix, "DATA/Harry/Coarfa SNP project/NHGRI Enrichment rev2/gwas_catalog_v1.0.2-associations_e100_r2020-10-07.tsv")), "\t", escape_double = FALSE, trim_ws = TRUE)

# Some SNPs are annotated weirdly and don't have entries in the position column
gwas_weird <- subset(gwas_cat, is.na(CHR_POS))
gwas_weird <- separate(gwas_weird, SNPS, into=(c("CHR_ID", "CHR_POS")), sep="[:_]")
gwas_weird <- gwas_weird[gwas_weird$CHR_ID %in% chrom_ids,]
gwas_weird$CHR_POS <- as.numeric(gwas_weird$CHR_POS)

# Give the GWAS IDs a good looking chromosome ID
gwas_cat$CHR_ID <- paste0("chr", gwas_cat$CHR_ID)

# Merge the weird ones back in
gwas_weird$CHR_POS <- as.character(gwas_weird$CHR_POS)
gwas_cat <- bind_rows(subset(gwas_cat, !is.na(CHR_ID)), gwas_weird)

# Get the mapped traits
gwas_cat <- separate(gwas_cat, MAPPED_TRAIT, into=c("trait_1", "trait_2", "trait_3"), sep=",")

# Apply a custom ID to each SNP based on its chromosome and position
gwas_cat$hm_index <- paste0("chr", gwas_cat$CHR_ID,"_", gwas_cat$CHR_POS)

# Clean up the SNPs more
gwas_sub <- subset(gwas_cat, select=c("CHR_ID", "CHR_POS", "SNPS", "STRONGEST SNP-RISK ALLELE","RISK ALLELE FREQUENCY", "P-VALUE", "hm_index", "trait_1", "trait_2", "trait_3"))
colnames(gwas_sub) <- c("chr", "pos", "snp", "risk_allele","raf", "p_val", "hm_index", "trait_1", "trait_2", "trait_3")

# Stack up SNPs with multiple traits so they can be merged
gwas_sub_1 <- subset(gwas_sub, select=c("chr", "pos", "snp", "risk_allele","raf", "p_val", "hm_index", "trait_1"))
colnames(gwas_sub_1)[8] <- "hm_trait"
gwas_sub_2 <- subset(gwas_sub, select=c("chr", "pos", "snp", "risk_allele","raf", "p_val", "hm_index", "trait_2"))
colnames(gwas_sub_2)[8] <- "hm_trait"
gwas_sub_3 <- subset(gwas_sub, select=c("chr", "pos", "snp", "risk_allele","raf", "p_val", "hm_index", "trait_3"))
colnames(gwas_sub_3)[8] <- "hm_trait"

gwas_sub_2 <- subset(gwas_sub_2, !is.na(hm_trait))
gwas_sub_3 <- subset(gwas_sub_3, !is.na(hm_trait))

# Stack up all the SNP-trait matches
gwas_sub <- bind_rows(gwas_sub_1, gwas_sub_2, gwas_sub_3)

# I've hand annotated all these trait categories with a scheme similar to Bonder et al., 2015
gwas_catalog_annotations <- fread((paste0(path_prefix, "DATA/Harry/Coarfa SNP project/NHGRI Enrichment rev2/hm custom gwas annotations rev2.csv")))
colnames(gwas_catalog_annotations)[3] <- "hm_trait"
hm_cat_ls <- unique(gwas_catalog_annotations$hm_category)
gwas_sub <- inner_join(gwas_sub, subset(gwas_catalog_annotations, select=c("hm_trait", "hm_category")), by="hm_trait")

# Sometimes the risk allele is actually the major allele (raf>0.5), but to keep things consistent we want to look at SNPs in terms of their minor allele freq
# So for any SNP with risk allele freq >0.5, compute the maf by subtracting it from 1
gwas_sub$raf <- as.numeric(gwas_sub$raf)
gwas_sub$risk_allele_status <- "A2"
gwas_sub$risk_allele_status[gwas_sub$raf>0.5] <- "A1"
gwas_sub$maf[gwas_sub$risk_allele_status=="A1"] <- 1-gwas_sub$raf[gwas_sub$raf>0.5]
gwas_sub$maf[gwas_sub$risk_allele_status=="A2"] <- gwas_sub$raf[gwas_sub$raf<0.5]

# We only want snps with more than 5% maf to be involved
gwas_sub <- subset(gwas_sub, maf>=0.05)

# Rename columns in the background list
colnames(freq_df)[1] <- "snp_id_38"
freq_df$start <- freq_df$pos
freq_df$end <- freq_df$pos

# Oh and no sex chromosomes
freq_df <- subset(freq_df, chr!="chrX" & chr!="chrY")
gwas_sub <- subset(gwas_sub, chr!="chrX" & chr!="chrY")

# Load all the Simes SNPs and use these to get Corsiv coordinates
simes_snps_all <- data.frame()


# Rob brought up trying this while disregarding tissue-specificity. So I'll just smush all the separate tissues together, then run the rest of the script
for (i in 1:length(tissue_list)){
  
  # Load the CoRSIV SNPs
  corsiv_snp <- fread(paste0(path_prefix, "DATA/Harry/Coarfa SNP project/NHGRI Enrichment rev2/simes_snps/indexed_distance_beta.", tissue_list[i], ".sig_fdr.0.05.csv"))
  
  # Separate the CoRSIV SNP's address
  corsiv_snp %>% separate(snp, c("chr", "pos")) -> corsiv_snp
  corsiv_snp$hm_index <- paste0(corsiv_snp$chr,"_", corsiv_snp$pos)
  corsiv_snp$pos <- as.numeric(corsiv_snp$pos)
  
  # Add to list of all simes snps
  simes_snps_all <- bind_rows(simes_snps_all, corsiv_snp)
  
}

simes_snps_all %>% subset(select=c("Chrom", "ChromStart", "ChromStop")) %>% distinct() -> corsiv_bound
colnames(corsiv_bound) <- c("chr", "start", "end")

# Load the corsivs and establish boundaries based on corsiv midpoint
corsiv_bound <- mutate(corsiv_bound, length=end-start, mid=start+0.5(length),  start_new=(mid-corsiv_window), end_new=mid+corsiv_window)
corsiv_bound <- subset(corsiv_bound, select=c("chr", "start_new", "end_new"))
colnames(corsiv_bound) <- c("chr", "start", "end")
corsiv_bound$start[which(corsiv_bound$start<0)] <- 0

# Some of these corsivs have boundaries that extend over the end of the chromosome
corsiv_bound <- left_join(corsiv_bound, hg38_chrom_lengths, by="chr")
corsiv_bound$end[corsiv_bound$end>corsiv_bound$chr_length] <- corsiv_bound$chr_length[corsiv_bound$end>corsiv_bound$chr_length]
corsiv_bound <- subset(corsiv_bound, start<end)


# Pick out all the background SNPs that are within the CORSIV boundaries
corsiv_bound_bg <- fuzzyjoin::genome_inner_join(corsiv_bound, subset(freq_df, select=c("chr", "start", "end", "pos", "snp")), by=c("chr", "start", "end"))
corsiv_bound_bg <- distinct(corsiv_bound_bg)
colnames(corsiv_bound_bg)[5:7] <- c("chr", "start", "end")

# Write this to disc so I don't need to regenerate this each time (it takes forever!)
fwrite(gwas_sub, "output_csv/nhgri gwas snps processed.csv")
fwrite(gwas_cat, "output_csv/nhgri catalogue processed.csv")
fwrite(corsiv_bound_bg, "output_csv/gtex snps flanking corsivs.csv")
