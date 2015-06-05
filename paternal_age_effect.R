#################################################################
# This file tests for associations between aneuploidy incidence #
# and paternal age.										                          #
#################################################################

library(ggplot2)
library(gridExtra)
library(Matching)
library(ppcor)
library(RCurl)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

data_filtered <- filterData(data)

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_blastomere <- callPloidy(data_blastomere)

####################

aneuploid_binom <- aneuploidyByCase(data_blastomere) 
paternal_age <- data.frame(aggregate(data_filtered$paternal_age ~ data_filtered$case, FUN = mean))
names(paternal_age) <- c("i", "paternal_age")
aneuploid_binom <- merge(aneuploid_binom, paternal_age, "i")

###########################
### Partial correlation ###
###########################

pcor.test((aneuploid_binom$aneuploid_1/(aneuploid_binom$aneuploid_0+aneuploid_binom$aneuploid_1)), aneuploid_binom$paternal_age, aneuploid_binom$maternal_age, method="spearman")
pcor.test((aneuploid_binom$aneuploid_1/(aneuploid_binom$aneuploid_0+aneuploid_binom$aneuploid_1)), aneuploid_binom$maternal_age, aneuploid_binom$paternal_age, method="spearman")

####################
# caliper = 0.1    #
####################

aneuploid_binom$father <- NA
aneuploid_binom[aneuploid_binom$paternal_age < median(aneuploid_binom$paternal_age),]$father <- 0
aneuploid_binom[aneuploid_binom$paternal_age >= median(aneuploid_binom$paternal_age),]$father <- 1

set.seed(1)

matches <- Match(Tr = aneuploid_binom$father, X = aneuploid_binom$maternal_age, replace = FALSE, ties = FALSE, caliper = 0.1)
aneuploid_binom_matched <- cbind(aneuploid_binom[matches$index.treated,], aneuploid_binom[matches$index.control,])

sum(aneuploid_binom_matched[,2]); sum(aneuploid_binom_matched[,3])
sum(aneuploid_binom_matched[,8]); sum(aneuploid_binom_matched[,9])

sum(aneuploid_binom_matched[,3])/sum(aneuploid_binom_matched[,2]+aneuploid_binom_matched[,3])
sum(aneuploid_binom_matched[,9])/sum(aneuploid_binom_matched[,8]+aneuploid_binom_matched[,9])

fisher.test(rbind(cbind(sum(aneuploid_binom_matched[,2]), sum(aneuploid_binom_matched[,3])), cbind(sum(aneuploid_binom_matched[,8]), sum(aneuploid_binom_matched[,9]))))

####################
# caliper = 0.01    #
####################

aneuploid_binom$father <- NA
aneuploid_binom[aneuploid_binom$paternal_age < median(aneuploid_binom$paternal_age),]$father <- 0
aneuploid_binom[aneuploid_binom$paternal_age >= median(aneuploid_binom$paternal_age),]$father <- 1

set.seed(1)

matches <- Match(Tr = aneuploid_binom$father, X = aneuploid_binom$maternal_age, replace = FALSE, ties = FALSE, caliper = 0.01)
aneuploid_binom_matched <- cbind(aneuploid_binom[matches$index.treated,], aneuploid_binom[matches$index.control,])

sum(aneuploid_binom_matched[,2]); sum(aneuploid_binom_matched[,3])
sum(aneuploid_binom_matched[,8]); sum(aneuploid_binom_matched[,9])

sum(aneuploid_binom_matched[,3])/sum(aneuploid_binom_matched[,2]+aneuploid_binom_matched[,3])
sum(aneuploid_binom_matched[,9])/sum(aneuploid_binom_matched[,8]+aneuploid_binom_matched[,9])

fisher.test(rbind(cbind(sum(aneuploid_binom_matched[,2]), sum(aneuploid_binom_matched[,3])), cbind(sum(aneuploid_binom_matched[,8]), sum(aneuploid_binom_matched[,9]))))

###########################
### GLM F-Test   	    ###
###########################

g0 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age, family = quasibinomial(link = "logit"))
g1 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2), family = quasibinomial(link = "logit"))

anova(g0, g1, test = "F")

g2 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3), family = quasibinomial(link = "logit"))

anova(g1, g2, test = "F")

g3 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3) + I(maternal_age^4), family = quasibinomial(link = "logit"))

anova(g2, g3, test = "F")

g4 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3) + paternal_age, family = quasibinomial(link = "logit"))

anova(g2, g4, test = "F")

###########################

donor <- data.frame(aggregate(data_filtered$egg_donor ~ data_filtered$case, FUN = mean))
names(donor) <- c("i", "donor")
aneuploid_binom <- merge(aneuploid_binom, donor, "i")

f1 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ paternal_age, family = quasibinomial(link = "logit"))
f2 <- glm(data = aneuploid_binom[aneuploid_binom$donor == 1,], formula = cbind(aneuploid_1, aneuploid_0) ~ paternal_age, family = quasibinomial(link = "logit"))

#########################################

data_te <- selectSampleType(data_filtered, TE)

data_te <- callPloidy(data_te)

####################

aneuploid_binom <- aneuploidyByCase(data_te) 
paternal_age <- data.frame(aggregate(data_filtered$paternal_age ~ data_filtered$case, FUN = mean))
names(paternal_age) <- c("i", "paternal_age")
aneuploid_binom <- merge(aneuploid_binom, paternal_age, "i")

###########################
### Partial correlation ###
###########################

pcor.test((aneuploid_binom$aneuploid_1/(aneuploid_binom$aneuploid_0+aneuploid_binom$aneuploid_1)), aneuploid_binom$paternal_age, aneuploid_binom$maternal_age, method="spearman")
pcor.test((aneuploid_binom$aneuploid_1/(aneuploid_binom$aneuploid_0+aneuploid_binom$aneuploid_1)), aneuploid_binom$maternal_age, aneuploid_binom$paternal_age, method="spearman")

####################
# caliper = 0.1    #
####################

aneuploid_binom$father <- NA
aneuploid_binom[aneuploid_binom$paternal_age < median(aneuploid_binom$paternal_age),]$father <- 0
aneuploid_binom[aneuploid_binom$paternal_age >= median(aneuploid_binom$paternal_age),]$father <- 1

set.seed(1)

matches <- Match(Tr = aneuploid_binom$father, X = aneuploid_binom$maternal_age, replace = FALSE, ties = FALSE, caliper = 0.1)
aneuploid_binom_matched <- cbind(aneuploid_binom[matches$index.treated,], aneuploid_binom[matches$index.control,])

sum(aneuploid_binom_matched[,2]); sum(aneuploid_binom_matched[,3])
sum(aneuploid_binom_matched[,8]); sum(aneuploid_binom_matched[,9])

sum(aneuploid_binom_matched[,3])/sum(aneuploid_binom_matched[,2]+aneuploid_binom_matched[,3])
sum(aneuploid_binom_matched[,9])/sum(aneuploid_binom_matched[,8]+aneuploid_binom_matched[,9])

fisher.test(rbind(cbind(sum(aneuploid_binom_matched[,2]), sum(aneuploid_binom_matched[,3])), cbind(sum(aneuploid_binom_matched[,8]), sum(aneuploid_binom_matched[,9]))))

###########################
### GLM F-Test   	    ###
###########################

g0 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age, family = quasibinomial(link = "logit"))
g1 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2), family = quasibinomial(link = "logit"))

anova(g0, g1, test = "F")

g2 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3), family = quasibinomial(link = "logit"))

anova(g1, g2, test = "F")

g3 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2) + paternal_age, family = quasibinomial(link = "logit"))

anova(g1, g3, test = "F")

###########################

donor <- data.frame(aggregate(data_filtered$egg_donor ~ data_filtered$case, FUN = mean))
names(donor) <- c("i", "donor")
aneuploid_binom <- merge(aneuploid_binom, donor, "i")

f1 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ paternal_age, family = quasibinomial(link = "logit"))
f2 <- glm(data = aneuploid_binom[aneuploid_binom$donor == 1,], formula = cbind(aneuploid_1, aneuploid_0) ~ paternal_age, family = quasibinomial(link = "logit"))
