####################################################################
# This file generates a plot of biopsies submitted by maternal age #
####################################################################

library(ggplot2)
library(gridExtra)
source("~/Desktop/aneuploidy_functions.R")

data <- read.table("~/Desktop/upload.csv", sep=",", header=T) # import the data

data_blastomere <- selectSampleType(data, blastomere)

data_te <- selectSampleType(data, TE)

blastomere_counts <- data.frame(table(data_blastomere$case))
names(blastomere_counts) <- c("case", "count")

te_counts <- data.frame(table(data_te$case))
names(te_counts) <- c("case", "count")

maternal_age <- data.frame(cbind(data[!duplicated(data$maternal_age),]$case, data[!duplicated(data$maternal_age),]$maternal_age))
names(maternal_age) <- c("case", "maternal_age")

blastomere_counts <- merge(maternal_age, blastomere_counts, "case")
te_counts <- merge(maternal_age, te_counts, "case")


samplesByAge <- function(data, label) {
  results <- data.frame(matrix(ncol = 5))
  minAge <- min(data$maternal_age[!is.na(data$maternal_age)])
  maxAge <- max(data$maternal_age[!is.na(data$maternal_age)])
  
  for (k in round(minAge):round(maxAge)) {	
    age_subset <- data[round(data$maternal_age) == k & !is.na(round(data$maternal_age)),]
    mean_samples <- mean(age_subset$count)
    se <- std(age_subset$count)
    samples <- c(mean_samples, se, k, label, nrow(age_subset))
    results <- rbind(results, samples)
  }
  results$X1 <- as.numeric(results$X1)
  results$X2 <- as.numeric(results$X2)
  results$X3 <- as.numeric(results$X3)
  results$X5 <- as.numeric(results$X5)
  return(results[-1,])
}

blastomere_count_results <- samplesByAge(blastomere_counts, "Day-3 blastomere")
te_count_results <- samplesByAge(te_counts, "Day-5 TE")

results_all_chroms <- rbind(blastomere_count_results, te_count_results)

limits<-aes(ymax=(X1+X2), ymin=(X1-X2))

ggplot(data = results_all_chroms, aes(x = X3, y = X1, col = factor(X4))) + 
  geom_errorbar(limits, width = 0.5) + geom_line() + geom_point() + xlab("Maternal Age") + ylab("Biopsies Submitted for Testing") + 
  theme(legend.position = "right") + theme_bw() + scale_color_discrete(name = "") + ylim(0, 17)
