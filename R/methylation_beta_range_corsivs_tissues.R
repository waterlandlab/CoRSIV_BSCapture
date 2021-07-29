library(ggpubr)
blood <- read.csv("../mQTL_results/blood_GTEx_mQTL.csv")
brain <- read.csv("../mQTL_results/brain_GTEx_mQTL.csv")
lung <- read.csv("../mQTL_results/lung_GTEx_mQTL.csv")
nerge <- read.csv("../mQTL_results/nerve_GTEx_mQTL.csv")
skin <- read.csv("../mQTL_results/skin_GTEx_mQTL.csv")
thyroid <- read.csv("../mQTL_results/thyroid_GTEx_mQTL.csv")

gghistogram(bins = 20,data,x ="methy_range",fill = "gray",color="black",xlab="Range in % Methylation",title = "% Methylation Range in CoRSIVs-Thyroid")+theme_classic(base_size = 25)

range(median(blood$methy_range),
median(brain$methy_range),
median(lung$methy_range),
median(skin$methy_range),
median(thyroid$methy_range))
