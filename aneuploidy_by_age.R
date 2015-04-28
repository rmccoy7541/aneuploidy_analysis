data_te <- selectSampleType(data_filtered, TE)

data_te <- callPloidy(data_TE)

#################################################################

results_all_chroms <- aneuploidyByAge(data_te, "Data")

aneuploid_binom <- aneuploidyByCase(data_te) # get summary of aneuploid / euploid embryos per case

g0 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age, family = quasibinomial(link = "logit"))

g1 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2), family = quasibinomial(link = "logit"))

anova(g0, g1, test="F")

g2 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3), family = quasibinomial(link = "logit"))

anova(g1, g2, test="F")

# Fit model to data
g3 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3) + I(maternal_age^4), family = quasibinomial(link = "logit"))

anova(g2, g3, test = "F")

# Use model to predict aneuploidy based on age
prediction_frame <- data.frame("maternal_age" = seq(20,48,1))
predicted <- data.frame(NA, NA, 20:48, "Fitted Values", NA)
names(predicted) <- c("X1", "X2", "X3", "X4", "X5")
predicted_values <- predict.glm(g1, newdata = prediction_frame, type = "response", se.fit = TRUE)
predicted$X1 <- predicted_values$fit
predicted$X2 <- predicted_values$se.fit
results_all_chroms <- rbind(results_all_chroms, predicted)

# Plot the results
results_all_chroms$X1 <- as.numeric(results_all_chroms$X1)
results_all_chroms$X2 <- as.numeric(results_all_chroms$X2)
results_all_chroms$X3 <- as.numeric(results_all_chroms$X3)
results_all_chroms$X5 <- as.numeric(results_all_chroms$X5)
limits <- aes(ymax = (X1 + X2), ymin = (X1 - X2))

ggplot(data = results_all_chroms[results_all_chroms$X5 > 9 | is.na(results_all_chroms$X5),], aes(x = X3, y = X1, col = factor(X4))) + geom_errorbar(limits, width = 0.5) + geom_line() + geom_point()  + xlab("Maternal Age") + ylab("Proportion Aneuploid Biopsies") + theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "")

# Test the effect of egg donor status
donor <- data.frame(aggregate(data_te$egg_donor ~ data_te$case, FUN = mean))
names(donor) <- c("i", "donor")
aneuploid_binom <- merge(aneuploid_binom, donor, "i")

g4 <- glm(data = aneuploid_binom, formula = cbind(aneuploid_1, aneuploid_0) ~ maternal_age + I(maternal_age^2) + I(maternal_age^3) + donor, family = quasibinomial(link = "logit"))

anova(g3, g4, test = "F")


results_te <- results_all_chroms
results_te$X6 <- "Day-3 TE"

##########################################################

results_all_chroms <- rbind(results_blastomere, results_te)

c <- ggplot(data = results_all_chroms[results_all_chroms$X5 > 9 | is.na(results_all_chroms$X5),], aes(x = X3, y = X1, col = factor(X4))) + geom_errorbar(limits, width = 0.5) + geom_line() + geom_point() + xlab("Maternal Age") + ylab("Proportion Aneuploid") + theme(legend.position = "right") + theme_bw() + scale_color_discrete(name = "")

c <- c + facet_grid(. ~ X6)
