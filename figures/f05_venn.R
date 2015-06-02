#################################################################
# This file generates venn diagram of complex aneuploidy        #
#################################################################

library(ggplot2)
library(gridExtra)
library(VennDiagram)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

data <- read.table("~/Desktop/aneuploidy_analysis-master/upload.csv", sep = ",", header = T) # import the data

data_filtered <- filterData(data)

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_blastomere <- callPloidy(data_blastomere)

maternal_trisomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_blastomere)))
for (i in 7:29) {
	new<-(data_blastomere[,i]=="H210" & data_blastomere[,i+69]==0)
	maternal_trisomy_frame[,i-6]<-new
}
paternal_trisomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_blastomere)))
for (i in 7:29) {
	new<-(data_blastomere[,i]=="H120" & data_blastomere[,i+69]==0)
	paternal_trisomy_frame[,i-6]<-new
}
maternal_monosomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_blastomere)))
for (i in 7:29) {
	new<-(data_blastomere[,i]=="H010" & data_blastomere[,i+69]==0)
	maternal_monosomy_frame[,i-6]<-new
}
paternal_monosomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_blastomere)))
for (i in 7:29) {
	new<-(data_blastomere[,i]=="H100" & data_blastomere[,i+69]==0)
	paternal_monosomy_frame[,i-6]<-new
}
nullisomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_blastomere)))
for (i in 7:29) {
	new<-(data_blastomere[,i]=="H000" & data_blastomere[,i+69]==0)
	nullisomy_frame[,i-6]<-new
}

area1 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
area2 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
area3 <- sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
area4 <- sum(apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
area5 <- sum(apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n12 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n13 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n14 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n15 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n23 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n24 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n25 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n34 <- sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n35 <- sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n45 <- sum(apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n123 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n124 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n125 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n134 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n135 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n145 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n234 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n235 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n245 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n345 <- sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n1234 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n1235 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n1245 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n1345 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n2345 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n12345 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)

countsA <- draw.quintuple.venn(area1 = area1, area2 = area2, area3 = area3, area4 = area4, area5 = area5, 
n12 = n12, n13 = n13, n14 = n14, n15 = n15, n23 = n23, n24 = n24, n25 = n25, n34 = n34, n35 = n35, 
n45 = n45, n123 = n123, n124 = n124, n125 = n125, n134 = n134, n135 = n135, n145 = n145, n234 = n234, 
n235 = n235, n245 = n245, n345 = n345, n1234 = n1234, n1235 = n1235, n1245 = n1245, n1345 = n1345, 
n2345 = n2345, n12345 = n12345, category=c("Maternal Trisomy", "Paternal Trisomy", "Maternal Monosomy", 
"Paternal Monosomy", "Nullisomy"), lwd=c(0,0,0,0,0), fill=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", 
"#ff7f00"),  fontfamily = "sans")


#################################################################

data_te <- selectSampleType(data_filtered, TE)

data_te <- callPloidy(data_te)


maternal_trisomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_te)))
for (i in 7:29) {
	new<-(data_te[,i]=="H210" & data_te[,i+69]==0)
	maternal_trisomy_frame[,i-6]<-new
}
paternal_trisomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_te)))
for (i in 7:29) {
	new<-(data_te[,i]=="H120" & data_te[,i+69]==0)
	paternal_trisomy_frame[,i-6]<-new
}
maternal_monosomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_te)))
for (i in 7:29) {
	new<-(data_te[,i]=="H010" & data_te[,i+69]==0)
	maternal_monosomy_frame[,i-6]<-new
}
paternal_monosomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_te)))
for (i in 7:29) {
	new<-(data_te[,i]=="H100" & data_te[,i+69]==0)
	paternal_monosomy_frame[,i-6]<-new
}
nullisomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_te)))
for (i in 7:29) {
	new<-(data_te[,i]=="H000" & data_te[,i+69]==0)
	nullisomy_frame[,i-6]<-new
}

area1 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
area2 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
area3 <- sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
area4 <- sum(apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
area5 <- sum(apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n12 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n13 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n14 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n15 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n23 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n24 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n25 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n34 <- sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n35 <- sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n45 <- sum(apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n123 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n124 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n125 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n134 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n135 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n145 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n234 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n235 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n245 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n345 <- sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n1234 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n1235 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n1245 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n1345 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n2345 <- sum(apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
n12345 <- sum(apply(maternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_trisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)

countsB <- draw.quintuple.venn(area1 = area1, area2 = area2, area3 = area3, area4 = area4, area5 = area5, 
n12 = n12, n13 = n13, n14 = n14, n15 = n15, n23 = n23, n24 = n24, n25 = n25, n34 = n34, n35 = n35, 
n45 = n45, n123 = n123, n124 = n124, n125 = n125, n134 = n134, n135 = n135, n145 = n145, n234 = n234, 
n235 = n235, n245 = n245, n345 = n345, n1234 = n1234, n1235 = n1235, n1245 = n1245, n1345 = n1345, 
n2345 = n2345, n12345 = n12345, category=c("Maternal Trisomy", "Paternal Trisomy", "Maternal Monosomy", 
"Paternal Monosomy", "Nullisomy"), lwd=c(0,0,0,0,0), fill=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", 
"#ff7f00"),  fontfamily = "sans")
