###NES distribution
library(ggplot2)
library(ggpubr)
load("/Users/XZX/thesis/addupFullNew3.RData")
load("/Users/XZX/thesis/addupCountNew3.RData")
distribution <- function(celltype,  method = "t.test", groupby = "Group", labely = 0.75){
  selected <- addupFullNew3[addupFullNew3$CellType == celltype,]
  p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
  p1 <- p1 + geom_boxplot() +
    facet_wrap(~ Group+Signature) +
    stat_compare_means(method = method,label = "p.signif",label.x = 1.5, label.y = labely) +
    labs(title = celltype, x = "Cluster Type", y = "NES")
  return(p1)
  print(p1)
}
categrz1 <- function(en_raw){
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
alldata <- function(nes, cluster_zhca_i){
  cluster <- cluster_zhca_i
  cluster$Var2 <- rownames(cluster)
  nes <- categrz1(nes)
  nes <- melt(nes)
  NES <- left_join(nes, cluster, by = c("variable" = "Var2"))
  return(NES)
}
load("/Users/XZX/thesis/Data_NES_immune_cells/G1_NES.RData")
load("/Users/XZX/thesis/Data_ClusterType/G1.RData")
NES1 <- alldata(nes, cluster_zhca_i)
NES1$group <- "G1"
load("/Users/XZX/thesis/Data_NES_immune_cells/G2_NES.RData")
load("/Users/XZX/thesis/Data_ClusterType/G2.RData")
NES2 <- alldata(nes, cluster_zhca_i)
NES2$group <- "G2"
load("/Users/XZX/thesis/Data_NES_immune_cells/G3_NES.RData")
load("/Users/XZX/thesis/Data_ClusterType/G3.RData")
NES3 <- alldata(nes, cluster_zhca_i)
NES3$group <- "G3"
load("/Users/XZX/thesis/Data_NES_immune_cells/G4_NES.RData")
load("/Users/XZX/thesis/Data_ClusterType/G4.RData")
NES4 <- alldata(nes, cluster_zhca_i)
NES4$group <- "G4"
NES <- rbind(NES1, NES2, NES3, NES4)
p1 <- distribution("Preadipocytes")

celltype = "Preadipocytes"
selected <- addupFullNew3[addupFullNew3$CellType == celltype,]
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_wrap(~Group+Signature) +
  stat_compare_means(method = "anova",label = "p.signif",label.x = 1.5, label.y = 0.75)+
  labs(title = celltype, x = "Cluster Type", y = "NES")

celltype = "Skeletal muscle"
selected <- NES[NES$CellType == celltype,]
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_grid(cols = vars(group), rows = vars(Signature), as.table = TRUE ) +
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5, label.y = 0.6) +
  labs(title = celltype, x = "Cluster Type", y = "NES")

Skeletal_muscle_ENCODE_1___457
Skeletal_muscle_FANTOM_2___461
selected <- NES[NES$Signature %in% c("mv_Endothelial_cells_FANTOM_1___355","mv_Endothelial_cells_FANTOM_2___356","mv_Endothelial_cells_FANTOM_3___357")]
selected <- filter(NES, Signature %in% c("Skeletal_muscle_ENCODE_1___457","Skeletal_muscle_FANTOM_2___461"))
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_grid(cols = vars(group), rows = vars(Signature), as.table = TRUE ) +
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.4, label.y = 0.67) +
  labs(title = celltype, x = "Cluster Type", y = "NES")

celltype <- "Endothelial cells"
selected <- NES[NES$Signature %in% c("mv_Endothelial_cells_FANTOM_1___355","mv_Endothelial_cells_FANTOM_2___356","mv_Endothelial_cells_FANTOM_3___357")]
selected <- filter(NES, Signature %in% c("mv_Endothelial_cells_FANTOM_1___355","mv_Endothelial_cells_FANTOM_3___357"))
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_grid(cols = vars(group), rows = vars(Signature), as.table = TRUE ) +
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.4, label.y = 0.67) +
  labs(title = celltype, x = "Cluster Type", y = "NES")

cell <- c("Preadipocytes", "MSC" ,"Fibroblasts", "Epithelial cells")
selected <- addupFullNew3[addupFullNew3$CellType =="CD8 Tcells",]
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_wrap(~Group) +
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5, label.y = 0.75)

selected <- addupFullNew2[addupFullNew2$CellType == "NKT ",]
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_wrap(~Group+TumorType) +
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5, label.y = 0.75)

cell <- "Macrophages"
selected <- addupFullNew3[addupFullNew3$CellType == cell,]
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_wrap(~Group) +
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5, label.y = 0.75) +
  labs(title = cell)

cell <- "NKT"
selected <- addupFullNew3[addupFullNew3$CellType == cell,]
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_wrap(~Group) +
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5, label.y = 0.75) +
  labs(title = cell)

cell <- "MSC"
selected <- addupFullNew3[addupFullNew3$CellType == cell,]
p1 <- ggplot(data = selected, mapping = aes(x = ClusterType, y = value))
p1 + geom_boxplot() +
  facet_wrap(~Group) +
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5, label.y = 0.65) +
  labs(title = cell)
