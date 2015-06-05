library(ggplot2)
library(gridExtra)
library(gtable)
library(RCurl)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

data_calls<-data[,7:29]; nrow(data_calls)
data_calls<-data_calls[(apply(data_calls, 1, function(x) sum(is.na(x)))<5),]; nrow(data_calls) # remove samples with 5+ no-calls
data_calls<-data_calls[(apply(data_calls, 1, function(x) sum(x[!is.na(x)]=="H000"))!=23),]; nrow(data_calls) # remove whole genome nullsomy
data_filtered<-data[row.names(data_calls),] #recover the original data frame, with all chrom. nullsomy and any chrom. no-calls filtered out

# create a new field indicating whole chromosome ploidy status
aneuploid_frame<-data.frame(matrix(ncol = 23, nrow= 42725))
for (i in 7:29) {
	new<-(data_filtered[,i]!="H110" & data_filtered[,i]!="H101" & data_filtered[,i+69]!=1)
	aneuploid_frame[,i-6]<-new
}
aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE))>0)
data_filtered$ploidy<-TRUE
data_filtered$ploidy[aneuploid_indicator]<-FALSE

####

chroms_affected<-apply(aneuploid_frame, 1, function(x) sum(x[!is.na(x)]==TRUE))
data_filtered$chroms_affected<-chroms_affected

blastomere<-hist(data_filtered[data_filtered$sample_type=="blastomere" & data_filtered$chroms_affected!=0,]$chroms_affected, breaks=seq(0,23,1))
te<-hist(data_filtered[data_filtered$sample_type=="TE" & data_filtered$chroms_affected!=0,]$chroms_affected, breaks=seq(0,23,1))

p<-te$counts/nrow(data[data$sample_type=="TE",])
q<-blastomere$counts/nrow(data[data$sample_type=="blastomere",])

results<-data.frame(rbind(cbind(seq(1,23,1), p, sqrt((p*(1-p))/nrow(data[data$sample_type=="TE",])), "Day-5 TE biopsies"), cbind(seq(1,23,1), q, sqrt((q*(1-q))/nrow(data[data$sample_type=="blastomere",])), "Day-3 blastomeres")))
names(results)<-c("chrom", "prop", "se", "sample_type")

results$chrom<-as.numeric(as.character(results$chrom))
results$prop<-as.numeric(as.character(results$prop))
results$se<-as.numeric(as.character(results$se))
limits<-aes(ymax=(prop+se), ymin=(prop-se))


a<-ggplot(data=results, aes(x=chrom, y=prop, fill=sample_type)) + geom_bar(stat="identity", position="dodge") + ylab("Prop. aneuploid samples") + xlab("Chromosomes affected") + theme(legend.position="right") + labs(fill="Sample type")# + geom_errorbar(limits, position=position_dodge(0.9), width=.5)

b<-ggplot(data=results, aes(x=chrom, y=prop, fill=sample_type)) + geom_bar(stat="identity", position="dodge") + ylab("Proportion Aneuploid Samples") + xlab("Number of Chromosomes Affected") + labs(fill="Sample type") + theme(legend.justification=c(0,0), legend.position=c(.25,0.55)) # + geom_errorbar(limits, position=position_dodge(0.9), width=.5) 


difference<-(-1*((results[results$sample_type=="Day-3 blastomeres",]$prop-results[results$sample_type=="Day-5 TE biopsies",]$prop)/results[results$sample_type=="Day-3 blastomeres",]$prop))

results2<-data.frame(cbind(1:23, difference))
names(results2)<-c("chrom", "difference")

c<-ggplot(data=results2, aes(x=chrom, y=difference, fill=factor(c(1, rep(2,22))))) + geom_bar(width=0.8, stat="identity", position="identity") + scale_fill_manual(values=c("#66c2a5", "indianred3")) + theme(legend.position="none") + ylab("Relative Difference") + xlab("Number of Chromosomes Affected")


a<-ggplot_gtable(ggplot_build(a))
b<-ggplot_gtable(ggplot_build(b))
c<-ggplot_gtable(ggplot_build(c))

c$widths<-b$widths


grid.newpage()
grid.arrange(b,c,nrow=2, heights=c(1,0.5))
