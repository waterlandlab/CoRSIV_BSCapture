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
library(openxlsx)

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


# I've hand annotated all these trait categories with a scheme similar to Bonder et al., 2015
gwas_catalog_annotations <- fread((paste0(path_prefix, "DATA/Harry/Coarfa SNP project/NHGRI Enrichment rev2/hm custom gwas annotations rev2.csv")))
colnames(gwas_catalog_annotations)[3] <- "hm_trait"
hm_cat_ls <- unique(gwas_catalog_annotations$hm_category)

# Load histogram data
all_hist <- data.frame()
all_stats <- data.frame()
for(i in 1:length(hm_cat_ls)){
  x <- fread(paste0("output_csv/background random match summary - ", hm_cat_ls[i], " - all tissues.csv"))
  x$trait_category <- hm_cat_ls[i]
  all_hist <- bind_rows(all_hist, x)
  
  y <- fread(paste0("output_csv/z test statistics - ", hm_cat_ls[i], ".csv"))
  all_stats <- bind_rows(all_stats, y)
}

all_stats <- mutate(all_stats, bootstrap_enrichment=n_corsiv_gwas_match/background_expected)
all_stats$p_adj <- p.adjust(all_stats$p_val, method="bonferroni")

# Plot a scatterplot
p <- ggplot(all_stats, aes(y=-log10(p_adj), x=bootstrap_enrichment, label=trait_category))+
  geom_hline(yintercept=-log10(0.05), linetype="dotted", colour="grey50")+
  geom_point(aes(fill=trait_category), shape=21, colour="black", size=3)+
  geom_text(hjust=0)+
  theme_hm()+
  scale_fill_brewer(palette="Set1")+
  theme(legend.position="none")+
  labs(y="-log10(p)", x="Enrichment")

func_ggsave(p, "point - log10p vs enrichment - all tissues", 4,4)

# Plot a load of histograms
p <- ggplot()+
  geom_histogram(data=all_hist, aes(x=overlap_sum, fill=trait_category), bins=35, alpha=0.75)+
  geom_point(data=all_stats, aes(x=n_corsiv_gwas_match, y=0), shape=23, fill="red4", size=3)+
  facet_wrap(~trait_category, nrow=2)+
  theme_hm()+
  theme(strip.background = element_blank())+
  lims(y=c(0,5000))+
  scale_fill_brewer(palette="Set1")+
  theme(legend.position="none")+
  labs(x="Matching GWAS SNPs", y="Count")

func_ggsave(p, "histogram - facetted enrichment - all tissues", 6,4)


out <- createWorkbook()

for (i in 1:length(hm_cat_ls)){
  temp_x <- fread(paste0("output_csv/gwas catalogue data on corsiv matched snps - 1000000 - ", hm_cat_ls[i], " - all tissues.csv"))
  addWorksheet(out, paste0("Table X - ", hm_cat_ls[i]))
  writeData(out, sheet=paste0("Table X - ", hm_cat_ls[i]), x=temp_x)
  
}

saveWorkbook(out, "supplemental table - corsiv overlapping nhgri snps.xlsx")
