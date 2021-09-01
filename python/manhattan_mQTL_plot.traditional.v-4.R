#!/usr/bin/env Rscript
library(ggplot2)
library(R.utils)

args = commandArgs(trailingOnly=TRUE)

input_file = args[1]
output_root= args[2]
pvalue_cutoff = as.double(args[3])
tissue = args[4]
# setwd("/Users/coarfa/Research/Waterland/Capture-RuiChen/manuscript/outlines/figure-3-support")
# input_file = "anno3.collated_annotated.xls"
# output_root = "anno3.collated_annotated"
# pvalue_cutoff = 0.001

printf("input %s output root %s pvalue cutoff %f tissue %s \n", input_file, output_root, pvalue_cutoff, tissue)

x_mQTL = read.table(input_file, header=T)

x_mQTL = x_mQTL[x_mQTL$pvalue<pvalue_cutoff,]
x_mQTL$minusLog10_pvalue = -log10(x_mQTL$pvalue)

x_mQTL$sign_beta = sign(x_mQTL$beta)
x_mQTL$beta_x_minusLog10_pvalue = x_mQTL$sign_beta * x_mQTL$minusLog10_pvalue

main_title = paste(tissue, x_mQTL[1,]$CoRSIV, x_mQTL[1,]$gene, sep=" ")

require(scales)
options(scipen=2000000)
p = ggplot(x_mQTL, aes(Distance, beta_x_minusLog10_pvalue)) + xlab("Distance from CoRSIV") + ylab("-log10(p-value)") + ggtitle(main_title)
p = p + geom_point(aes(colour = factor(sign_beta)), size = 0.5) + xlim(-1000000,1000000)  + theme_minimal() + expand_limits(x = 0, y = 0) + theme_classic() 
p = p+geom_hline(yintercept=0)
p = p + theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_rect(colour = "black", size=1))
p
# p+scale_x_continuous(labels = comma)
output_linear = paste(output_root,".gene.",x_mQTL[1,]$gene,".linear.pdf", sep="")
ggsave(output_linear, width=9, height=3, dpi=150)
#output_linear = paste(output_root,".linear.jpg", sep="")
#ggsave(output_linear, width=9, height=3, dpi=150)


# next do logarithmic scale
x_mQTL$Distance_add_1 = x_mQTL$Distance
x_mQTL$Distance_add_1[x_mQTL$Distance>=0] = x_mQTL$Distance[x_mQTL$Distance>=0] +1
x_mQTL$logDistance[x_mQTL$Distance>=0]=log2(x_mQTL$Distance_add_1[x_mQTL$Distance>=0])
x_mQTL$Distance_add_1[x_mQTL$Distance<0] = x_mQTL$Distance[x_mQTL$Distance<0] -1
x_mQTL$logDistance[x_mQTL$Distance<0]=-log2(abs(x_mQTL$Distance_add_1[x_mQTL$Distance<0]))

p = ggplot(x_mQTL, aes(logDistance, beta_x_minusLog10_pvalue)) + xlab("Distance from CoRSIV (log2)") + ylab("-log10(p-value)") + ggtitle(main_title)
p = p + geom_point(aes(colour = factor(sign_beta)), size = 0.5) + xlim(-log2(1000000), log2(1000000))  + theme_minimal() + expand_limits(x = 0, y = 0)+ theme_classic() 
p = p+geom_hline(yintercept=0)
p = p + theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_rect(colour = "black", size=1))
p
output_log2 = paste(output_root,".gene.",x_mQTL[1,]$gene,".log2.pdf", sep="")
ggsave(output_log2, width=9, height=3, dpi=150)
#output_log2 = paste(output_root,".log2.jpg", sep="")
#ggsave(output_log2, width=9, height=3, dpi=150)


