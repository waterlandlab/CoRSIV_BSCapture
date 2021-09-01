#!/usr/bin/env Rscript

printf <- function(...)print(sprintf(...))

library(MatrixEQTL)
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR

args = commandArgs(trailingOnly=TRUE)
snp_file = args[1]
methylation_file = args[2]
output_file = args[3]

cat("got command line arguments ",  str(args) , "\n")
printf("SNP file: %s methylation file %s output file %s\n",snp_file, methylation_file, output_file)


SNP_file_name = snp_file 
expression_file_name = methylation_file 
covariates_file_name=character()
# setwd("/Users/coarfa/Research/Waterland/Capture-RuiChen/manuscript/SNV-vs-Methylation/example")
output_file_name = output_file 
pvOutputThreshold = 1;
errorCovariance = numeric();

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
# cvrt$LoadFile(covariates_file_name);

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
