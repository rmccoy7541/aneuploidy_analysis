#################################################################
# This file tests rates of aneuploidy versus various referral   #
# reasons.                                                      #
#################################################################

library(ggplot2)
library(gridExtra)
library(RCurl)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

data_filtered <- filterData(data)

data_filtered <- data_filtered[data_filtered$egg_donor == 0,] # remove egg donor cases when testing for associations with referral reasons

data_blastomere <- selectSampleType(data_filtered, blastomere)
data_te <- selectSampleType(data_filtered, TE)

data_blastomere <- callPloidy(data_blastomere)
data_te <- callPloidy(data_te)
data_blastomere <- callMeiotic(data_blastomere)
data_te <- callMeiotic(data_te)
data_blastomere <- callMitotic(data_blastomere)
data_te <- callMitotic(data_te)

data_blastomere <- callMaternalTriploidy(data_blastomere)
data_te <- callMaternalTriploidy(data_te)

#################################################################

blastomere_aneuploidy_counts <- aneuploidyByCase(data_blastomere)
names(blastomere_aneuploidy_counts) <- c("case", "euploid", "aneuploid")

te_aneuploidy_counts <- aneuploidyByCase(data_te)
names(te_aneuploidy_counts) <- c("case", "euploid", "aneuploid")

#################################################################

maternal_age <- data.frame(cbind(data[!duplicated(data$case),]$case, data[!duplicated(data$case),]$maternal_age))
names(maternal_age) <- c("case", "maternal_age")

referral_reasons <- read.table("~/Desktop/aneuploidy_analysis-master/referral_reasons.csv", sep = ",", header = T)

referral_reasons <- merge(referral_reasons, maternal_age, "case")

#################################################################

referral_reasons_blastomere <- merge(referral_reasons, blastomere_aneuploidy_counts, "case")

g1 <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ maternal_age, family = quasibinomial(link = "logit"))
g2 <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2), family = quasibinomial(link = "logit"))

anova(g1, g2, test = "F")

g3 <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3), family = quasibinomial(link = "logit"))

anova(g2, g3, test = "F")

aneuploidy_blastomere <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + repeat_pregnancy_loss + previous_ivf_failure + male_factor + unexplained_infertility + translocation + previous_aneuploidy, family = quasibinomial(link = "logit"))

#################################################################

referral_reasons_te <- merge(referral_reasons, te_aneuploidy_counts, "case")

g1 <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age, family = quasibinomial(link = "logit"))
g2 <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2), family = quasibinomial(link = "logit"))

anova(g1, g2, test = "F")

g3 <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3), family = quasibinomial(link = "logit"))

anova(g2, g3, test = "F")

aneuploidy_te <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + repeat_pregnancy_loss + previous_ivf_failure + male_factor + unexplained_infertility + translocation + previous_aneuploidy, family = quasibinomial(link = "logit"))

#################################################################
# Test mitotic-origin aneuploidies separately                   #
#################################################################

blastomere_mitotic_counts <- mitoticByCase(data_blastomere)
names(blastomere_mitotic_counts) <- c("case", "euploid", "aneuploid")

te_mitotic_counts <- mitoticByCase(data_te)
names(te_mitotic_counts) <- c("case", "euploid", "aneuploid")

#################################################################

referral_reasons_blastomere <- merge(referral_reasons, blastomere_mitotic_counts, "case")

mitotic_blastomere <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ repeat_pregnancy_loss + previous_ivf_failure + male_factor + unexplained_infertility + translocation + previous_aneuploidy, family = quasibinomial(link = "logit"))

#################################################################

referral_reasons_te <- merge(referral_reasons, te_mitotic_counts, "case")

mitotic_te <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ repeat_pregnancy_loss + previous_ivf_failure + male_factor + unexplained_infertility + translocation + previous_aneuploidy, family = quasibinomial(link = "logit"))

#################################################################
# Test meiotic-origin aneuploidies separately                   #
#################################################################

blastomere_meiotic_counts <- meioticByCase(data_blastomere)
names(blastomere_meiotic_counts) <- c("case", "euploid", "aneuploid")

te_meiotic_counts <- meioticByCase(data_te)
names(te_meiotic_counts) <- c("case", "euploid", "aneuploid")

#################################################################

referral_reasons_blastomere <- merge(referral_reasons, blastomere_meiotic_counts, "case")

g1 <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ maternal_age, family = quasibinomial(link = "logit"))
g2 <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2), family = quasibinomial(link = "logit"))

anova(g1, g2, test = "F")

g3 <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3), family = quasibinomial(link = "logit"))

anova(g2, g3, test = "F")

meiotic_blastomere <- glm(data = referral_reasons_blastomere, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + repeat_pregnancy_loss + previous_ivf_failure + male_factor + unexplained_infertility + translocation + previous_aneuploidy, family = quasibinomial(link = "logit"))

#################################################################

referral_reasons_te <- merge(referral_reasons, te_meiotic_counts, "case")

g1 <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age, family = quasibinomial(link = "logit"))
g2 <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2), family = quasibinomial(link = "logit"))

anova(g1, g2, test = "F")

g3 <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3), family = quasibinomial(link = "logit"))

anova(g2, g3, test = "F")

g4 <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3) + I(maternal_age^4), family = quasibinomial(link = "logit"))

anova(g3, g4, test = "F")

g5 <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3) + I(maternal_age^4) + I(maternal_age^5), family = quasibinomial(link = "logit"))

anova(g4, g5, test = "F")

meiotic_te <- glm(data = referral_reasons_te, formula = cbind(euploid, aneuploid) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3) + I(maternal_age^4) + repeat_pregnancy_loss + previous_ivf_failure + male_factor + unexplained_infertility + translocation + previous_aneuploidy, family = quasibinomial(link = "logit"))

int_aneuploidy_blastomere <- confint(aneuploidy_blastomere)
int_aneuploidy_te <- confint(aneuploidy_te)
int_mitotic_blastomere <- confint(mitotic_blastomere)
int_mitotic_te <- confint(mitotic_te)
int_meiotic_blastomere <- confint(meiotic_blastomere)
int_meiotic_te <- confint(meiotic_te)

aneuploidy_blastomere <- summary(aneuploidy_blastomere)
aneuploidy_te <- summary(aneuploidy_te)
meiotic_blastomere <- summary(meiotic_blastomere)
meiotic_te <- summary(meiotic_te)
mitotic_blastomere <- summary(mitotic_blastomere)
mitotic_te <- summary(mitotic_te)

results_frame <- data.frame(rbind(
	cbind(attr(aneuploidy_blastomere$term, "term.labels")[3:8], aneuploidy_blastomere$coefficients[4:nrow(aneuploidy_blastomere$coefficients),1], aneuploidy_blastomere$coefficients[4:nrow(aneuploidy_blastomere$coefficients),2], "Aneuploidy", "Day-3 blastomeres"), 
	cbind(attr(aneuploidy_te$term, "term.labels")[3:8], aneuploidy_te$coefficients[4:nrow(aneuploidy_te$coefficients),1], aneuploidy_te$coefficients[4:nrow(aneuploidy_te$coefficients),2], "Aneuploidy", "Day-5 TE biopsies"),
	cbind(attr(meiotic_blastomere$term, "term.labels")[3:8], meiotic_blastomere$coefficients[4:nrow(meiotic_blastomere$coefficients),1], meiotic_blastomere$coefficients[4:nrow(meiotic_blastomere$coefficients),2], "Meiotic", "Day-3 blastomeres"), 
	cbind(attr(meiotic_te$term, "term.labels")[5:10], meiotic_te$coefficients[6:nrow(meiotic_te$coefficients),1], meiotic_te$coefficients[6:nrow(meiotic_te$coefficients),2], "Meiotic", "Day-5 TE biopsies"),
	cbind(attr(mitotic_blastomere$term, "term.labels"), mitotic_blastomere$coefficients[2:nrow(mitotic_blastomere$coefficients),1], mitotic_blastomere$coefficients[2:nrow(mitotic_blastomere$coefficients),2], "Mitotic", "Day-3 blastomeres"), 
	cbind(attr(mitotic_te$term, "term.labels"), mitotic_te$coefficients[2:nrow(mitotic_te$coefficients),1], mitotic_te$coefficients[2:nrow(mitotic_te$coefficients),2], "Mitotic", "Day-5 TE biopsies")
))

names(results_frame) <- c("terms", "beta", "se", "aneuploidy", "dev_stage")
results_frame$beta <- as.numeric(as.character(results_frame$beta))
results_frame$se <- as.numeric(as.character(results_frame$se))

limits <- aes(xmax = beta + se, xmin = beta - se)
p <- ggplot(data = results_frame, aes(x = beta, y = factor(terms), col = factor(terms))) + geom_point() + geom_errorbarh(limits, height = 0)
p + facet_grid(aneuploidy ~ dev_stage) + ylab('Referral Reason') + xlab('Coefficient') + geom_vline(xintercept = 0, linetype = "dotted") + theme(legend.position="none")

results_frame$beta <- (-1) * results_frame$beta
results_frame$or <- exp(1) ^ (results_frame$beta)

results_frame <- cbind(results_frame, rbind(int_aneuploidy_blastomere[4:nrow(int_aneuploidy_blastomere),], int_aneuploidy_te[4:nrow(int_aneuploidy_te),], int_meiotic_blastomere[4:nrow(int_meiotic_blastomere),], int_meiotic_te[6:nrow(int_meiotic_te),], int_mitotic_blastomere[2:nrow(int_mitotic_blastomere),], int_mitotic_te[2:nrow(int_mitotic_te),]))

results_frame$ul <- exp(1) ^ (results_frame[,7] * (-1))
results_frame$ll <- exp(1) ^ (results_frame[,8] * (-1))

limits <- aes(xmax = ul, xmin = ll)
p <- ggplot(data = results_frame, aes(x = or, y = factor(terms), col = factor(terms))) + geom_point() + geom_errorbarh(limits, height = 0)
p + facet_grid(aneuploidy ~ dev_stage) + ylab('Referral Reason') + xlab('Odds Ratio, 95% CI') + geom_vline(xintercept = 1, linetype = "dotted") + theme(legend.position="none")


q <- ggplot(data = results_frame[(results_frame$terms == "translocation" | results_frame$terms == "previous_ivf_failure" | results_frame$terms == "repeat_pregnancy_loss"),], aes(x = or, y = factor(terms), col = factor(terms))) + theme_bw() + geom_point() + geom_errorbarh(limits, height = 0)
q + facet_grid(aneuploidy ~ dev_stage) + ylab('Referral Reason') + xlab('Odds Ratio, 95% CI') + geom_vline(xintercept = 1, linetype = "dotted") + theme(legend.position="none") 
