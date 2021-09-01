#!/usr/bin/env Rscript
library(rsq)
library(modEvA)

printf <- function(...)print(sprintf(...))

args = commandArgs(trailingOnly=TRUE)
corsiv_support_file = args[1]
glm_output = args[2]

printf("Input file %s glm_output %s \n", corsiv_support_file, glm_output)
x_corsiv = read.table(corsiv_support_file, header=T, row.names=1)

xx = as.numeric(x_corsiv[1,])
yy = as.numeric(x_corsiv[3,])
my_glm = glm(yy~xx)
p_value = coef(summary(my_glm))[,4][[2]]
rsq_r_value = rsq(my_glm)
dsquared_r_value = Dsquared(my_glm)
glm_df = data.frame(p_value, rsq_r_value, dsquared_r_value)
write.csv(glm_df, file="my1stglm.txt")
