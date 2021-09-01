#!/usr/bin/env Rscript
library(rsq)
library(modEvA)

printf <- function(...)print(sprintf(...))

args = commandArgs(trailingOnly=TRUE)
corsiv_support_file = args[1]
glm_output = args[2]

printf("Input file %s glm_output %s \n", corsiv_support_file, glm_output)
x_corsiv = read.table(corsiv_support_file, header=T, row.names=1)

my_glm_0 = glm(x_corsiv$methylation ~ x_corsiv$allele_sum)
rsq_r_value_0 = rsq(my_glm_0)
dsquared_r_value_0 = Dsquared(my_glm_0)
ct_0 = cor.test(my_glm_0$y, my_glm_0$fitted.values)
pvalue_0 = ct_0$p.value

my_glm_3 = glm(x_corsiv$methylation ~ x_corsiv$allele1_sum + x_corsiv$allele2_sum + x_corsiv$allele_sum)
rsq_r_value_3 = rsq(my_glm_3)
dsquared_r_value_3 = Dsquared(my_glm_3)
ct_3 = cor.test(my_glm_3$y, my_glm_3$fitted.values)
pvalue_3= ct_3$p.value

my_glm_1 = glm(x_corsiv$methylation ~ x_corsiv$allele1_index + x_corsiv$allele2_index)
rsq_r_value_1 = rsq(my_glm_1)
dsquared_r_value_1 = Dsquared(my_glm_1)
ct_1 = cor.test(my_glm_1$y, my_glm_1$fitted.values)
pvalue_1 = ct_1$p.value

my_glm_2 = glm(x_corsiv$methylation ~ x_corsiv$allele1_index + x_corsiv$allele2_index + x_corsiv$allele1_sum + x_corsiv$allele2_sum + x_corsiv$allele_sum)
rsq_r_value_2 = rsq(my_glm_2)
dsquared_r_value_2 = Dsquared(my_glm_2)
ct_2 = cor.test(my_glm_1$y, my_glm_1$fitted.values)
pvalue_2 = ct_2$p.value

glm_df = data.frame(rsq_r_value_0, pvalue_0, 
                    rsq_r_value_1, pvalue_1, 
                    rsq_r_value_2, pvalue_2,
                    rsq_r_value_3, pvalue_3
                    )
write.csv(glm_df, file=glm_output)
