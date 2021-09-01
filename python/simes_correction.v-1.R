#!/usr/bin/env Rscript

printf <- function(...)print(sprintf(...))

#################################################################
## Script to perform the Simes procedure (Simes, 1986)
#################################################################

simes = function (p) {
        o = order(p, decreasing = TRUE)
        ro = order(o)
        lp = length(p)
        pmin(1, (p[o] * lp/(lp:1))[ro])
}


#################################################################
## Example
#################################################################

# set.seed(1)
# PVals = runif(10)
# simes(PVals)
## In the two-step approach, the minimum is taken as the p-value for a CpG
# min(simes(PVals))

args = commandArgs(trailingOnly=TRUE)
mQTL_file = args[1]
mQTL_w_simes = args[2]

printf("mQTL %s output %s\n",mQTL_file, mQTL_w_simes)

x = read.csv(mQTL_file)

head(x)

x$q_simes = simes(x$pValue)

write.csv(x, file=mQTL_w_simes)


