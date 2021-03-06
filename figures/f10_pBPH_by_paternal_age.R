#################################################################
# This file generates a plot of proportion of pBPH aneuploidies # 
# versus paternal age.                                          #
#################################################################

library(ggplot2)
library(gridExtra)
library(RCurl)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

data_filtered <- filterData(data)

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_blastomere <- callPloidy(data_blastomere)

data_blastomere <- callPaternalMeiotic(data_blastomere)

data_blastomere <- callMitotic(data_blastomere)

#################################################################

results_all_chroms <- paternalMeioticByAge(data_blastomere, "Data")


# Plot the results
results_all_chroms$X1 <- as.numeric(results_all_chroms$X1)
results_all_chroms$X2 <- as.numeric(results_all_chroms$X2)
results_all_chroms$X3 <- as.numeric(results_all_chroms$X3)
results_all_chroms$X5 <- as.numeric(results_all_chroms$X5)
limits <- aes(ymax = (X1 + X2), ymin = (X1 - X2))

results_blastomere <- results_all_chroms
results_blastomere$X6 <- "Day-3 Blastomere"

########################################

data_te <- selectSampleType(data_filtered, TE)

data_te <- callPloidy(data_te)

data_te <- callPaternalMeiotic(data_te)

data_te <- callMitotic(data_te)

#################################################################

results_all_chroms <- paternalMeioticByAge(data_te, "Data")

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

c <- ggplot(data = results_all_chroms[results_all_chroms$X5 > 9 | is.na(results_all_chroms$X5),], aes(x = X3, y = X1)) + 
geom_errorbar(limits, width = 0.5) + geom_line() + geom_point() + xlab("Paternal Age") + 
ylab("Prop. with Paternal BPH Error") + theme(legend.position = "right") + theme_bw() + scale_color_discrete(name = "")

c <- c + facet_grid(. ~ X6)

##########################################################

meiotic_blastomere <- meioticByCase(data_blastomere)
meiotic_te <- meioticByCase(data_te)

paternal_age <- data.frame(aggregate(data$paternal_age ~ data$case, FUN = mean))
names(paternal_age) <- c("case", "paternal_age")
meiotic_blastomere <- merge(meiotic_blastomere, paternal_age, "case")
meiotic_te <- merge(meiotic_te, paternal_age, "case")

summary(glm(data = meiotic_blastomere, formula = cbind(euploid, aneuploid) ~ paternal_age, family = quasibinomial(link = "logit")))
summary(glm(data = meiotic_te, formula = cbind(euploid, aneuploid) ~ paternal_age, family = quasibinomial(link = "logit")))
