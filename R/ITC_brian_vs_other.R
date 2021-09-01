library(dplyr)
library(tidyr)
library(ggpubr)

corsiv_methY_data <- read.table(file = "../CoRSIV Capture methylation_data/CoRSIV_Meth_data.txt",header = T, sep = "\t")
corsiv_methY_data <- corsiv_methY_data[
  with(corsiv_methY_data, order(CHR, STR,decreasing = F)),
]
grep("Brain",names(corsiv_methY_data))
grep("Blood",names(corsiv_methY_data))




pdf("blood_vs_other_scatterPlots_orderd.pdf",paper = "special",height = 25,width = 25 )
par(mfrow=c(5,5))
for(i in 1:2340){

  X <- data.frame(t(corsiv_methY_data[i,grep("Blood",names(corsiv_methY_data))]))
  X$sample_name <- as.character(rownames(X))
  X <- X %>% separate(sample_name, c("A","B","C"),"__")
  
  for(tissue in c("Brain","Thyroid","Skin","Lung","Tibial")){
    #tissue="Brain"
    Y <- data.frame(t(corsiv_methY_data[i,grep(tissue,names(corsiv_methY_data))]))
    Y$sample_name <- as.character(rownames(Y))
    Y <- Y %>% separate(sample_name, c("A","B","C"),"__")
    
    merged_df <- merge(X,Y,by="B")
    par(mar=c(5,5,5,5))
    mod1 = lm(merged_df[,5]~merged_df[,2], data = merged_df)
    modsum = summary(mod1)
    r2 = modsum$adj.r.squared
    mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
    
    CoRSIV_coord <- paste(strsplit(corsiv_methY_data[i,1],split = ";")[[1]][1],strsplit(corsiv_methY_data[i,1],split = ";")[[1]][2],strsplit(corsiv_methY_data[i,1],split = ";")[[1]][3],sep="_")
    plot(col.lab="red", cex.lab=2,cex.axis = 2,frame.plot = T,    #  for the xlab and ylab
         col="black",merged_df[,2],merged_df[,5],pch=16,xlab = "Blood % Methylation",ylab = paste0(tissue," % Methylation"),main=CoRSIV_coord,xlim = c(0,1),ylim = c(0,1),cex = 1.5)
    text(x = 0.2, y = 1, labels = mylabel)
    my.p = modsum$coefficients[2,4]
    text(x = 0.2, y = 0.9, labels = bquote(italic(P) == .(format(my.p, digits = 3))))
    
  }
    
}
dev.off()   
    
  
  # ##
  # Y <- data.frame(t(corsiv_methY_data[i,grep("Thyroid",names(corsiv_methY_data))]))
  # Y$sample_name <- as.character(rownames(Y))
  # Y <- Y %>% separate(sample_name, c("A","B","C"),"__")
  # 
  # merged_df <- merge(X,Y,by="B")
  # plot(col.lab="red", cex.lab=2,cex.axis = 2,   #  for the xlab and ylab
  #      col="black",merged_df[,2],merged_df[,5],pch=16,xlab = "Brain % Methylation",ylab = "Thyroid  % Methylation",xlim = c(0,1),ylim = c(0,1),cex = 1.5)
  # 
  # ##
  # Y <- data.frame(t(corsiv_methY_data[i,grep("Skin",names(corsiv_methY_data))]))
  # Y$sample_name <- as.character(rownames(Y))
  # Y <- Y %>% separate(sample_name, c("A","B","C"),"__")
  # 
  # merged_df <- merge(X,Y,by="B")
  # plot(col.lab="red", cex.lab=2,cex.axis = 2,   #  for the xlab and ylab
  #      col="black",merged_df[,2],merged_df[,5],pch=16,xlab = "Brain % Methylation",ylab = "Skin % Methylation",xlim = c(0,1),ylim = c(0,1),cex = 1.5)
  # 
  # ###
  # Y <- data.frame(t(corsiv_methY_data[i,grep("Lung",names(corsiv_methY_data))]))
  # Y$sample_name <- as.character(rownames(Y))
  # Y <- Y %>% separate(sample_name, c("A","B","C"),"__")
  # 
  # merged_df <- merge(X,Y,by="B")
  # plot(col.lab="red", cex.lab=2,cex.axis = 2,    #  for the xlab and ylab
  #      col="black",merged_df[,2],merged_df[,5],pch=16,xlab = "Brain % Methylation",ylab = "Lung % Methylation",xlim = c(0,1),ylim = c(0,1),cex = 1.5)
  # 
  # ##
  # Y <- data.frame(t(corsiv_methY_data[i,grep("Tibial",names(corsiv_methY_data))]))
  # Y$sample_name <- as.character(rownames(Y))
  # Y <- Y %>% separate(sample_name, c("A","B","C"),"__")
  # 
  # merged_df <- merge(X,Y,by="B")
  # plot(col.lab="red", cex.lab=2,cex.axis = 2,    #  for the xlab and ylab
  #      col="black",merged_df[,2],merged_df[,5],pch=16,xlab = "Brain % Methylation",ylab = "Tibial % Methylation",xlim = c(0,1),ylim = c(0,1),cex = 1.5)
  

