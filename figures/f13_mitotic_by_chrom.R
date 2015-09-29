##################################################################
# This file generates a plot of per-chromosome rates of putative #
# mitotic-origin aneuploidies.                                   #
##################################################################

library(ggplot2)
library(gridExtra)
library(RCurl)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

data_filtered <- filterData(data)

data_filtered <- callPloidy(data_filtered)
data_filtered <- callMaternalMeiotic(data_filtered)
data_filtered <- callMitotic(data_filtered)

#################################################################

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_blastomere <- callPloidy(data_blastomere)

data_blastomere <- callMitotic(data_blastomere)

blastomere_results <- mitoticByChromosome(data_blastomere) # calculate chrom-specific aneuploidy rates

blastomere_results$type <- "Blastomere, Day 3"

#################################################################

data_te <- selectSampleType(data_filtered, TE)

data_te <- callPloidy(data_te)

data_te <- callMitotic(data_te)

te_results <- mitoticByChromosome(data_te) # calculate chrom-specific aneuploidy rates

te_results$type  <- "Trophectoderm, Day 5"

#################################################################

results <- rbind(blastomere_results, te_results) # combine sample type results

results$type <- factor(results$type, levels=c("Blastomere, Day 3", "Trophectoderm, Day 5", "Miscarriage"))
limits <- aes(ymax = p + se, ymin = p - se)

a <- ggplot(data = results[order(results$chrom),], aes(y = p, x = chrom, fill = factor(type))) + 
geom_bar(stat = "identity", width = 0.8, position = "dodge") + 
geom_errorbar(limits, width = 0.25, position = position_dodge(0.8)) + 
scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", 
"17", "18", "19", "20", "21", "22", "X/Y")) + xlab('Chromosome') + ylab('Prop. with Mitotic Error') + 
theme(axis.text.x = element_text(angle = 65, hjust = 1, size=8))  + theme_bw() + theme(legend.position="top") + 
scale_fill_discrete("")

# plot correlations

limits <- aes(xmax = p + se, xmin = p - se, ymin = results[results$type == "Trophectoderm, Day 5",]$p, 
ymax = results[results$type == "Trophectoderm, Day 5",]$p)
limits2 <- aes(ymax = p + se, ymin = p - se, xmin=results[results$type == "Blastomere, Day 3",]$p, 
xmax=results[results$type == "Blastomere, Day 3",]$p)

b <- ggplot(data = results, aes(x = results[results$type == "Blastomere, Day 3",]$p, 
y = results[results$type == "Trophectoderm, Day 5",]$p)) + theme_bw() + 
geom_errorbarh(limits, data = results[results$type == "Blastomere, Day 3",]) + 
geom_errorbar(limits2, data = results[results$type == "Trophectoderm, Day 5",]) +  
geom_text(color = "gray", aes(label = as.character(c(1:22, as.character("X/Y"))))) + 
xlab("Blastomere, Day 3") + ylab("Trophectoderm, Day 5")

cor.test(results[results$type == "Blastomere, Day 3",]$p, results[results$type == "Trophectoderm, Day 5",]$p)

c <- ggplot(data = data_filtered[data_filtered$meiotic == FALSE & data_filtered$mitotic == TRUE & 
data_filtered$chroms_affected > 0,], aes(x = chroms_affected, fill = factor(sample_type))) + theme_bw() + 
geom_histogram(binwidth = 1, position = "dodge") + xlab('Total Affected Chromosomes') + ylab('Number of Samples') + 
theme(legend.position="none")

grid.arrange(a, arrangeGrob(b, c, nrow = 1, widths = c(0.51, 0.5)), nrow = 2) 
