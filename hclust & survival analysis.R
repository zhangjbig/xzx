#-------for hclust and survival analysis------#
path<-"/Users/XZX/thesis"
library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(pheatmap)
library(ggplotify)
categrz <- function(en_raw){
  dataframe <- as.data.frame(en_raw)
  dataframe$sgntr <- row.names(dataframe)
  library(stringr)
  for (i in 1:nrow(dataframe)){
    signtr_i <- str_split(row.names(dataframe)[i], "_\\d")[[1]]
    type_i <- str_c(signtr_i[1], sep = "_")
    dataframe$type[i] <- type_i
  }
  for (i in c("HPCA",   "IRIS",   "ENCODE",   "FANTOM",   "NOVERSHTERN",   "BLUEPRINT")){
    dataframe$type <- gsub(i,"", dataframe$type, fixed=TRUE)
  }
  dataframe$type <- gsub("_","", dataframe$type)
  dataframe
}
prepare <- function(nes, signcell){
  sign <- subset(signcell, signcell[,2] <= 0.05)
  nes <- nes[rownames(nes) %in% rownames(sign),]
  nest <- t(nes)
  zhc <- scale(nest, center = TRUE, scale = TRUE)
  return(zhc)
}


##-------input data-------##
load("/Users/XZX/thesis/Data_clinic_info/G4_clinicinfo.RData")##########clinicinfo
load("/Users/XZX/thesis/Data_NES_immune_cells_matched/G4_NES.RData")
load("/Users/XZX/thesis/Data_coxreg_immune cells/G4_cox_reg.RData")
Group <- "Group 4"
zhc <- prepare(nes, en_raw)
#title <- "Group 1 Re (Group 5)"
#title <- "Group 4"
#sign <- categrz(t(zhc))

##--------cluster--------##
i <- 2
d <- dist(zhc, method = "euclidean")
zhca <- hclust(d, method = "complete")
ClusterType <- cutree(zhca, k=2)
cluster_zhca_i <- as.data.frame(ClusterType)
cluster_zhca_i$ClusterType <- gsub(1,"Long",cluster_zhca_i$ClusterType)
cluster_zhca_i$ClusterType <- gsub(2,"Short",cluster_zhca_i$ClusterType)
#write.csv(cluster_zhca_i, file = "/Users/XZX/thesis/G1_Pri.csv")

#------heatmap------#
#pdf(file="myplot1.pdf", height=15, width=40)
colors=list(ClusterType=c("Long"="#1B9E77", "Short"="#D95F02"))
p <- pheatmap(t(zhc), border_color = NA,cellwidth = 5,cellheight = 5,treeheight_row = 0, main = Group, 
              annotation_col = cluster_zhca_i, annotation_colors = colors, fontsize_row = 5, 
              show_colnames = FALSE,fontsize_col = 5)
p4 <- as.ggplot(p)

#dev.off()

#------KM------#
SurvAnlys <-  merge(clinicinfo, cluster_zhca_i, by=0, all=FALSE)
fit <- survfit(Surv(time, status) ~ ClusterType, data = SurvAnlys)
diff <- pairwise_survdiff(Surv(time, status) ~ ClusterType, data = SurvAnlys)
km4 <- ggsurvplot(fit, title = Group,pval = TRUE, palette = c("#1B9E77", "#D95F02"), 
                 legend.labs=c("Long","Short"), legend.title="Cluster Type", 
                 xlab = "Survival Time (months)")
                 #risk.table = TRUE, tables.height = 0.2,tables.theme = theme_cleantable()
km4$plot
fit

library(cowplot)
plot_grid(km1$plot,km2$plot,km3$plot,km4$plot,ncol = 2,nrow=2)
plot_grid(p3, p4, nrow=2)
library(customLayout)

ggarrange(clusthm, KM,labels = c("A", "B"),ncol = 2, nrow = 1)

#col.names = c("Samples","CellSignatures","CellTypes","Long","Short","pValue")
#total <- data.frame(row.names = c("Group 1 Pri","Group 1 Re","Group 2","Group 3","Group 4","Group 6 Pri", "Group 6 Re"))

#total[1,] <- list(nrow(zhc), ncol(zhc), "CellTypes", fit[["strata"]][["ClusterType=Long"]],
                  #fit[["strata"]][["ClusterType=Short"]], diff[["p.value"]])
group <- list()
group$Samples<- nrow(zhc)
group$CellSignatures <- ncol(zhc)
group$CellTypes <- length(unique(sign[,ncol(sign)]))
group$Long <-fit[["strata"]][["ClusterType=Long"]]
group$Short <- fit[["strata"]][["ClusterType=Short"]]
group$pValue <- diff[["p.value"]]
group <- as.data.frame(group)
rownames(group) <- title
ifelse(total == '', total<-group, total<-rbind(total,group))
  
colnames(total) = c("Samples","CellSignatures","CellTypes","Long","Short","pValue")  
#save(total, file = "/Users/XZX/thesis/total.RData")
readr::write_csv(total, path="total.csv")

library(corrplot)
res_cor <- cor(total)
corrplot(corr=res_cor)


