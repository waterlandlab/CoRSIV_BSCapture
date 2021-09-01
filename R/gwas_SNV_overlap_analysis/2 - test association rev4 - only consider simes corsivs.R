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

corsiv_bound_bg <- fread("output_csv/gtex snps flanking corsivs.csv")
corsiv_bound_bg <- mutate(corsiv_bound_bg, corsiv_id=paste0(chr.x,":",start.x,"-", end.x)) # This corsiv ID doesn't match Chathura's since it involves the extended windows, but it doesn't matter since we're just using it to keep things organized

gwas_sub <- fread("output_csv/nhgri gwas snps processed.csv")
gwas_sub$chr <- stringr::str_replace(gwas_sub$chr, "chrchr", "chr")

corsiv_bound_bg <- subset(corsiv_bound_bg, select=c("chr", "pos", "snp", "corsiv_id"))

gwas_cat <- fread("output_csv/nhgri catalogue processed.csv")


result_df <- data.frame()
simes_snps_all <- data.frame()


# Rob brought up trying this while disregarding tissue-specificity. 
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

simes_snps_all <- distinct(simes_snps_all, hm_index, .keep_all=TRUE)
corsiv_snp <- simes_snps_all # The script expects the Simes snps to be in this dataframe

hm_cat_ls <- c("various", "metabolic", "neurological", "anthropometric", "immune", "cardiovascular", "cancer", "hematological")
corsiv_ls <- unique(corsiv_bound_bg$corsiv_id)

for (k in 1:length(hm_cat_ls)){
  
  # Pull out the GWAS hits that match the current annotation category
  gwas_sub_c <- subset(gwas_sub, hm_category==hm_cat_ls[k])
  gwas_sub_c <- distinct(gwas_sub_c, snp, .keep_all=TRUE)
  gwas_sub_c$pos <- as.numeric(gwas_sub_c$pos)
  gwas_sub_c$gwas_overlap <- 1
  
  # Now for each corsiv, we want to pick 1000 random SNPs from within its boundaries
  # For testing:
  #corsiv_ls_sub <- sample(corsiv_ls, 1000, replace=FALSE) # 
  #corsiv_bound_bg_sub <- corsiv_bound_bg[corsiv_id %in% corsiv_ls_sub,]
  
  # This is the heart of the operation, Dplyr speeds this up massively. For each corsiv, grab N SNVs (with replacement)
  # Then label each SNV with a number matching the 'iteration' (the first SNP grabbed for each corsiv is considered the first iteration, and so forth)
  system.time(corsiv_bound_bg %>% group_by(corsiv_id) %>% sample_n(10000, replace=TRUE) %>% mutate(iteration=seq(from=1, to=10000, by=1)) -> random_snps)
  random_snps <- left_join(random_snps, subset(gwas_sub_c, select=c("chr", "pos", "gwas_overlap")), by=c("chr", "pos")) # See how many overlap trait-associated SNVs
  random_snps$gwas_overlap[is.na(random_snps$gwas_overlap)] <- 0 
  random_snps %>% group_by(iteration) %>% dplyr::summarize(overlap_sum=sum(gwas_overlap)) -> background_stat
  #hist(background_stat$overlap_sum)
  
  # Now let's check overlap with our gwas snps real quick
  corsiv_gwas_match <- inner_join(gwas_sub_c, corsiv_snp, by=c("chr", "pos"))
  corsiv_gwas_match$hm_index <- paste0(corsiv_gwas_match$chr,"_", corsiv_gwas_match$pos)
  
  # Now we can fold the matching SNPs back into original GWAS data
  match_gwas_snp <- left_join(corsiv_gwas_match, gwas_cat, by="hm_index")
  match_gwas_snp$trait_1 <- as.character(match_gwas_snp$trait_1)
  match_gwas_snp$trait_2 <- as.character(match_gwas_snp$trait_2)
  match_gwas_snp$trait_3 <- as.character(match_gwas_snp$trait_3)
  
  match_gwas_snp <- subset(match_gwas_snp, select=c("chr", "pos", "snp", "risk_allele", "risk_allele_status", "hm_category", "hm_trait"))
  
  fwrite(match_gwas_snp, paste0("output_csv/gwas catalogue data on corsiv matched snps - 1000000 - ", hm_cat_ls[k], " - all tissues", ".csv"))
  
  # Let's save the null histogram for later too
  fwrite(background_stat, paste0("output_csv/background random match summary - ", hm_cat_ls[k], " - all tissues.csv"))
  
  fwrite(background_random, paste0("output_csv/background random matches - ", hm_cat_ls[k], " - all tissues.csv"))
  
  # Now do a Z-test on the corsiv matches
  background_random_matches_z <- scale(background_stat$overlap_sum)
  background_expected <- mean(background_stat$overlap_sum)
  
  # Z-transform the actual overlap on the same scale as the null distribution
  z_stat <- scale(nrow(corsiv_gwas_match), center=background_expected, scale=attr(background_random_matches_z,"scaled:scale")) 
  Z <- distributions3::Normal(0, 1)
  z_test <- 1 - distributions3::cdf(Z, abs(z_stat)) + distributions3::cdf(Z, -abs(z_stat))
  
  # Add the basic counts to a growing dataframe
  temp_df <- data.frame(hm_cat_ls[k], nrow(corsiv_gwas_match), nrow(corsiv_snp), background_expected, z_test[1,])
  colnames(temp_df) <- c("trait_category","n_corsiv_gwas_match", "n_corsiv_snp",  "background_expected", "p_val")
  
  # Add to results
  fwrite(temp_df, paste0("output_csv/z test statistics - ", hm_cat_ls[k], ".csv"))
         
}


