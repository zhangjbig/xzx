#--------------------------#
#-------NES score----------#


path<-"/Users/XZX/Desktop/thesis/Data"

#input parameters and data
library(readxl) 
library(plyr)
#input cell markers (each row is gene signature for a cell type)
markers <- read_excel(paste(path,"/curatedMarkers_20191113.xlsx",sep=""), col_names = FALSE)

#input expression profiling (rows are genes,columns are samples)
geneExp_pri <- read_excel(paste(path,"/geneExp.pri.xlsx",sep=""),na = "NaN")
gene<-gsub('^.|.$', '', geneExp_pri$...1)  #match the character at the start of the string (^.) or the end of the string (.$) and replace it with ''
exp<-""
exp<-data.frame(geneExp_pri[,2:dim(geneExp_pri)[2]])  #nrow0 <- function(x) dim(x)[1],ncol0 <- function(x) dim(x)[2]
rownames(exp)<-gene

#input samples and delete null samples in exp
sample<-gsub('^.|.$', '', colnames(geneExp_pri))
colnames(exp)<-sample[2:length(sample)]
exp<-exp[,!is.na(exp[1,])]  #delete null samples in exp

#run NES method
nes<-NULL
for (i in 1:nrow(markers)){
  en_raw<-''
  markers_i<-markers[i,!is.na(as.character(markers[i,]))]   #markers_i is a new data frame deleting null genes in cell_i
  cell_i<-as.character(markers_i[1])   #get the cell name of markers_i as char
  markers_i<-as.character(markers_i[2:length(markers_i)])  #get the marker genes of cell_i as char


  markers_exp<-''
  markers_exp<-intersect(rownames(exp),markers_i)  #the matched genes in a cell type
  outside_markers_exp<-setdiff(rownames(exp),markers_i)  #the unmatched gene in a cell type

  if (length(markers_exp)>=3){
    for (exp_i in 1 : dim(exp)[2]){  #traverse every sample in exp, every sample is retrieved as exp_i
      tmp<-''
      tmp <- wilcox.test(exp[markers_exp,exp_i], exp[outside_markers_exp,exp_i], 
                         alternative = "greater")  #compare the significance in difference btw gene expression of a sample in immune cells(or not)
      ans<-''
      ans$statistic <- tmp$statistic; names(ans$statistic) <- NULL  #store the value into ans and remove the name
      ans$nes <- tmp$statistic/length(markers_exp)/length(outside_markers_exp) - 
        (length(markers_exp) + 1)/(2*length(outside_markers_exp)); names(ans$nes) <- NULL
      tmp_nes<-''
      tmp_nes<-matrix(ans$nes,nrow=1)  #why need to mark "nrow=1"
      colnames(tmp_nes)<-colnames(exp)[exp_i]
      if (en_raw == ''){
        en_raw<-tmp_nes
      }else{
        en_raw<-cbind(en_raw,tmp_nes)
      }
    }
    rownames(en_raw)<-cell_i
    nes<-rbind(nes,en_raw)
  }

}

#to save data, please change the dirctory
save(nes,file= paste(path,"/NES_geneExp_priXZX.RData",sep=""))

