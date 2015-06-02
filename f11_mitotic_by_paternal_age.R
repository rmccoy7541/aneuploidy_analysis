#################################################################
# This file generates a plot of the proportion of putative      #
# mitotic-origin aneuploidies versus paternal age.              #
#################################################################

library(ggplot2)
library(gridExtra)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

data <- read.table("~/Desktop/aneuploidy_analysis-master/upload.csv", sep=",", header=T) # import the data

data_filtered <- filterData(data)

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_blastomere <- callPloidy(data_blastomere)

data_blastomere <- callPaternalMeiotic(data_blastomere)

data_blastomere <- callMitotic(data_blastomere)


#################################################################

results_all_chroms <- mitoticByPaternalAge(data_blastomere, "Data")


# Plot the results
results_all_chroms$X1 <- as.numeric(results_all_chroms$X1)
results_all_chroms$X2 <- as.numeric(results_all_chroms$X2)
results_all_chroms$X3 <- as.numeric(results_all_chroms$X3)
results_all_chroms$X5 <- as.numeric(results_all_chroms$X5)
limits <- aes(ymax = (X1 + X2), ymin = (X1 - X2))

ggplot(data = results_all_chroms[results_all_chroms$X5 > 9 | is.na(results_all_chroms$X5),], aes(x = X3, y = X1)) + geom_errorbar(limits, width = 0.5) + geom_line() + geom_point()  + xlab("Paternal Age") + ylab("Prop. Blastomeres with BPH Aneuploidy") + theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "")


results_blastomere <- results_all_chroms
results_blastomere$X6 <- "Day-3 Blastomere"

########################################

data_te <- selectSampleType(data_filtered, TE)

data_te <- callPloidy(data_te)

data_te <- callPaternalMeiotic(data_te)

data_te <- callMitotic(data_te)


#################################################################

results_all_chroms <- mitoticByPaternalAge(data_te, "Data")


# Plot the results
results_all_chroms$X1 <- as.numeric(results_all_chroms$X1)
results_all_chroms$X2 <- as.numeric(results_all_chroms$X2)
results_all_chroms$X3 <- as.numeric(results_all_chroms$X3)
results_all_chroms$X5 <- as.numeric(results_all_chroms$X5)
limits <- aes(ymax = (X1 + X2), ymin = (X1 - X2))

results_te <- results_all_chroms
results_te$X6 <- "Day-5 TE"

##########################################################

results_all_chroms <- rbind(results_blastomere, results_te)

c <- ggplot(data = results_all_chroms[results_all_chroms$X5 > 9 | is.na(results_all_chroms$X5),], aes(x = X3, y = X1)) + geom_errorbar(limits, width = 0.5) + geom_line() + geom_point() + xlab("Paternal Age") + ylab("Prop. with Mitotic Aneuploidy") + theme(legend.position = "right") + theme_bw() + scale_color_discrete(name = "")

c <- c + facet_grid(. ~ X6)
