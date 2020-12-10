#------cox regression and graph-----#
library(survival)
library(survminer)
library(pheatmap)
organize <- function(nes, clinicinfo){
  valid <- intersect(rownames(clinicinfo), colnames(nes))
  clinicinfo <- clinicinfo[rownames(clinicinfo) %in% valid,]
  nes <- nes[colnames(nes) %in% valid,]
  para <- list(clinicinfo, nes)
  return(para)
}
coxreg <- function(nes, clinicinfo){
  nest <- t(nes)
  total <- cbind(nest, clinicinfo)
  en_raw <- data.frame()
  total <- as.data.frame(total)
  for (i in 1: ncol(nest)){
    cox_i <- coxph(Surv(time, status) ~ total[,i], data = total)
    en_raw[i,1] <- summary(cox_i)$coefficients[, 2]
    en_raw[i,2] <- summary(cox_i)$coefficients[, 5]
  }
  colnames(en_raw) <- c("Hazard Ratio","p_value")
  rownames(en_raw) <- colnames(nest)
  return(en_raw)
}

prepare <- function(nes, signcell){
  sign <- subset(signcell, signcell[,2] <= 0.05)
  nes <- nes[rownames(nes) %in% rownames(sign),]
  nest <- t(nes)
  zhc <- scale(nest, center = TRUE, scale = TRUE)
  return(zhc)
}

#---------groups---------#
load("/Users/XZX/thesis/Data_clinic_info/G1_clinicinfo.RData")##########clinicinfo
load("/Users/XZX/thesis/Data_NES_immune_cells_matched/G1_NES.RData")
valid <- intersect(rownames(clinicinfo), colnames(nes))
clinicinfo <- clinicinfo[rownames(clinicinfo) %in% valid,]
nes <- nes[,colnames(nes) %in% valid]
en_raw <- coxreg(nes, clinicinfo)
save(en_raw, file = "/Users/XZX/thesis/Data_coxreg_immune cells/new/G1_cox_reg.RData")
save(clinicinfo, file = "/Users/XZX/thesis/Data_clinic_info/G1_clinicinfo.RData")
save(nes, file = "/Users/XZX/thesis/Data_NES_immune_cells/G1_NES.RData")

load("/Users/XZX/thesis/Data_clinic_info/G2_clinicinfo.RData")##########clinicinfo
load("/Users/XZX/thesis/Data_NES_immune_cells_matched/G2_NES.RData")
valid <- intersect(rownames(clinicinfo), colnames(nes))
clinicinfo <- clinicinfo[rownames(clinicinfo) %in% valid,]
nes <- nes[,colnames(nes) %in% valid]
en_raw <- coxreg(nes, clinicinfo)
save(en_raw, file = "/Users/XZX/thesis/Data_coxreg_immune cells/new/G2_cox_reg.RData")
save(clinicinfo, file = "/Users/XZX/thesis/Data_clinic_info/G2_clinicinfo.RData")
save(nes, file = "/Users/XZX/thesis/Data_NES_immune_cells/G2_NES.RData")

load("/Users/XZX/thesis/Data_clinic_info/G3_clinicinfo.RData")##########clinicinfo
load("/Users/XZX/thesis/Data_NES_immune_cells_matched/G3_NES.RData")
clinicinfo <- clinicinfo[clinicinfo$CENSORING_STATUS!="NA",]
valid <- intersect(rownames(clinicinfo), colnames(nes))
clinicinfo <- clinicinfo[rownames(clinicinfo) %in% valid,]
nes <- nes[,colnames(nes) %in% valid]
en_raw <- coxreg(nes, clinicinfo)
save(en_raw, file = "/Users/XZX/thesis/Data_coxreg_immune cells/new/G3_cox_reg.RData")
save(clinicinfo, file = "/Users/XZX/thesis/Data_clinic_info/G3_clinicinfo.RData")
save(nes, file = "/Users/XZX/thesis/Data_NES_immune_cells/G3_NES.RData")

load("/Users/XZX/thesis/Data_clinic_info/G4_clinicinfo.RData")##########clinicinfo
load("/Users/XZX/thesis/Data_NES_immune_cells_matched/G4_NES.RData")
valid <- intersect(rownames(clinicinfo), colnames(nes))
clinicinfo <- clinicinfo[rownames(clinicinfo) %in% valid,]
nes <- nes[,colnames(nes) %in% valid]
en_raw <- coxreg(nes, clinicinfo)
save(en_raw, file = "/Users/XZX/thesis/Data_coxreg_immune cells/new/G4_cox_reg.RData")
save(clinicinfo, file = "/Users/XZX/thesis/Data_clinic_info/G4_clinicinfo.RData")
save(nes, file = "/Users/XZX/thesis/Data_NES_immune_cells/G4_NES.RData")

total <- ''
#---------graph---------#
load("/Users/XZX/thesis/Data_clinic_info/G4_clinicinfo.RData")##########clinicinfo
load("/Users/XZX/thesis/Data_NES_immune_cells_matched/G4_NES.RData")
load("/Users/XZX/thesis/Data_coxreg_immune cells/new/G4_cox_reg.RData")
load("/Users/XZX/thesis/Data_ClusterType/G4.RData")
Group <- "Group 4"
##--------cluster--------##
zhc <- prepare(nes, en_raw)
d <- dist(zhc, method = "euclidean")
zhca <- hclust(d, method = "complete")
#ClusterType <- cutree(zhca, k=2)
#cluster_zhca_i <- as.data.frame(ClusterType)
#cluster_zhca_i$ClusterType <- gsub(2,"Short",cluster_zhca_i$ClusterType)
#cluster_zhca_i$ClusterType <- gsub(1,"Long",cluster_zhca_i$ClusterType)
#save(cluster_zhca_i, file = "/Users/XZX/thesis/Data_ClusterType/G4.RData")

write.csv(cluster_zhca_i, file = "/Users/XZX/thesis/Data_ClusterType/G3.csv")
#------heatmap------#
#pdf(file="myplot1.pdf", height=15, width=40)
colors=list(ClusterType=c("Long"="#1B9E77", "Short"="#D95F02"))
p <- pheatmap(t(zhc), border_color = NA,
              cellwidth = 5,cellheight = 5,
              treeheight_row = 0,
              cutree_cols =2, 
              main = Group, 
              annotation_col = cluster_zhca_i, 
              annotation_colors = colors,
              fontsize_row = 5, 
              show_colnames = FALSE,fontsize_col = 5)
p4 <- as.ggplot(p)

#------KM------#

SurvAnlys <-  merge(clinicinfo, cluster_zhca_i, by=0, all=FALSE)
fit <- survfit(Surv(time, status) ~ ClusterType, data = SurvAnlys)
diff <- pairwise_survdiff(Surv(time, status) ~ ClusterType, data = SurvAnlys)
km4 <- ggsurvplot(fit, 
                  title = Group,
                  pval = TRUE, palette = c("#1B9E77", "#D95F02"), 
                  legend.labs=c("Long","Short"), legend.title="Cluster Type", 
                  xlab = "Survival Time (months)")
#risk.table = TRUE, tables.height = 0.2,tables.theme = theme_cleantable()
km4$plot
fit

library(cowplot)
plot_grid(km1$plot,km2$plot,km3$plot,km4$plot,ncol = 4,nrow=1)
#func categrz
celltype <- function(en_raw){
  dataframe <- as.data.frame(en_raw)
  dataframe$Signature <- row.names(dataframe)
  library(stringr)
  for (i in 1:nrow(dataframe)){
    signtr_i <- str_split(row.names(dataframe)[i], "_\\d")[[1]]
    type_i <- str_c(signtr_i[1], sep = "_")
    dataframe$CellType[i] <- type_i
  }
  for (i in c("HPCA",   "IRIS",   "ENCODE",   "FANTOM",   "NOVERSHTERN",   "BLUEPRINT")){
    dataframe$CellType <- gsub(i,"", dataframe$CellType, fixed=TRUE)
  }
  dataframe$CellType <- gsub("_+"," ", dataframe$CellType)
  dataframe$CellType <- str_trim(dataframe$CellType, side = "right")
  dataframe$name <- gsub("___\\d+","",dataframe$Signature)
  dataframe
  #for (i in 1:nrow(dataframe)){
  # id_i <- str_split(row.names(dataframe)[i], "_+")
  #if (j %in% dataframe$CellType){dataframe$id <- gsub(j,"", id_i, fixed=TRUE)}
  #dataframe$id[i] <- str_c(id_i, sep="_")
  #dataframe$id <- gsub("_+"," ", dataframe$CellType)}
}
hcn <- celltype(t(zhc))
group <- list()
group$Group <- Group
group$Samples <- ncol(hcn)
group$Signatures <- nrow(hcn)
group$`Cell Types` <- length(unique(hcn$CellType))
group$Long <- fit[["strata"]][["ClusterType=Long"]]
group$Short <- fit[["strata"]][["ClusterType=Short"]]
group$`p value` <- as.numeric(diff[["p.value"]])
group <- as.data.frame(group)
rownames(group) <- Group
ifelse(total == '', total<-group, total<-rbind(total,group))