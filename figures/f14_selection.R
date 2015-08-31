library(ggplot2)
library(gridExtra)
library(gtable)
library(RCurl)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

data_filtered <- filterData(data)

data_filtered <- callPloidy(data_filtered)

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_te <- selectSampleType(data_filtered, TE)


blastomere <- hist(data_blastomere[data_blastomere$chroms_affected != 0,]$chroms_affected, breaks=seq(0, 23, 1))
te <- hist(data_te[data_te$chroms_affected != 0,]$chroms_affected, breaks=seq(0, 23, 1))

p <- te$counts / nrow(data_te)
q <- blastomere$counts / nrow(data_blastomere)

results <- data.frame(rbind(cbind(seq(1,23,1), p, sqrt((p * (1 - p)) / nrow(data_te)), "Day-5 TE biopsies"), cbind(seq(1,23,1), q, sqrt((q * (1 - q)) / nrow(data_blastomere)), "Day-3 blastomeres")))
names(results)<-c("chrom", "prop", "se", "sample_type")

results$chrom <- as.numeric(as.character(results$chrom))
results$prop <- as.numeric(as.character(results$prop))
results$se <- as.numeric(as.character(results$se))
limits <- aes(ymax = (prop + se), ymin = (prop - se))


a <- ggplot(data=results, aes(x=chrom, y=prop, fill=sample_type)) + geom_bar(stat="identity", position="dodge") + ylab("Prop. aneuploid samples") + xlab("Chromosomes affected") + theme(legend.position="right") + labs(fill="Sample type")# + geom_errorbar(limits, position=position_dodge(0.9), width=.5)

b <- ggplot(data=results, aes(x=chrom, y=prop, fill=sample_type)) + geom_bar(stat="identity", position="dodge") + theme_bw() + ylab("Proportion Non-euploid Samples") + xlab("Number of Chromosomes Affected") + labs(fill="Sample type") + theme(legend.justification=c(0,0), legend.position=c(.25,0.55)) # + geom_errorbar(limits, position=position_dodge(0.9), width=.5) 


difference<-(-1*((results[results$sample_type=="Day-3 blastomeres",]$prop-results[results$sample_type=="Day-5 TE biopsies",]$prop)/results[results$sample_type=="Day-3 blastomeres",]$prop))

results2<-data.frame(cbind(1:23, difference))
names(results2)<-c("chrom", "difference")

c <- ggplot(data=results2, aes(x=chrom, y=difference, fill=factor(c(1, rep(2,22))))) + theme_bw() + geom_bar(width=0.8, stat="identity", position="identity") + scale_fill_manual(values=c("#66c2a5", "indianred3")) + theme(legend.position="none") + ylab("Relative Difference") + xlab("Number of Chromosomes Affected")


a<-ggplot_gtable(ggplot_build(a))
b<-ggplot_gtable(ggplot_build(b))
c<-ggplot_gtable(ggplot_build(c))

c$widths<-b$widths


grid.newpage()
grid.arrange(b,c,nrow=2, heights=c(1,0.5))
