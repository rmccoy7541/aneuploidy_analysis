#################################################################
# This file generates a plot of maternal vs. paternal and gain  #
# vs. loss aneuploidy versus the total number of chromosomes    #
# affected.                                                     #
#################################################################

library(ggplot2)
library(gridExtra)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

data_filtered <- filterData(data)

data_filtered <- callPloidy(data_filtered)

#################################

maternal_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data_filtered)))
for (i in 7:29) {
  new <- (data_filtered[,i] == "H200" | data_filtered[,i] == "H020" | data_filtered[,i] == "H010" | data_filtered[,i] == "H001" | data_filtered[,i] == "H000" | data_filtered[,i] == "H210" | data_filtered[,i] == "H201" | data_filtered[,i] == "H021") & data_filtered[,i + 69] != 1
  maternal_frame[,i - 6] <- new
}

paternal_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data_filtered)))
for (i in 7:29) {
  new <- (data_filtered[,i] == "H200" | data_filtered[,i] == "H020" | data_filtered[,i] == "H100" | data_filtered[,i] == "H000" | data_filtered[,i] == "H102" | data_filtered[,i] == "H120" | data_filtered[,i] == "H201" | data_filtered[,i] == "H021" | data_filtered[,i] == "H111") & data_filtered[,i + 69] != 1
  paternal_frame[,i - 6] <- new
}

maternal <- apply(maternal_frame, 1, function(x) sum(x[!is.na(x)] == TRUE))
data_filtered$maternal <- maternal

paternal <- apply(paternal_frame, 1, function(x) sum(x[!is.na(x)] == TRUE))
data_filtered$paternal <- paternal

data_filtered$prop_maternal_affected <- data_filtered$maternal / (data_filtered$maternal + data_filtered$paternal)

#################################

gains_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data_filtered)))
for (i in 7:29) {
  new <- (data_filtered[,i] == "H200" | data_filtered[,i] == "H020" | data_filtered[,i] == "H102" | data_filtered[,i] == "H210" | data_filtered[,i] == "H120" | data_filtered[,i] == "H201" | data_filtered[,i] == "H021" | data_filtered[,i] == "H111") & data_filtered[,i + 69] != 1
  gains_frame[,i - 6] <- new
}

gains <- apply(gains_frame, 1, function(x) sum(x[!is.na(x)] == TRUE))
data_filtered$gains <- gains

losses_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data_filtered)))
for (i in 7:29) {
  new <- (data_filtered[,i] == "H200" | data_filtered[,i] == "H020" | data_filtered[,i] == "H010" | data_filtered[,i] == "H100" | data_filtered[,i] == "H001" | data_filtered[,i] == "H000") & data_filtered[,i + 69] != 1
  losses_frame[,i - 6] <- new
}

losses <- apply(losses_frame, 1, function(x) sum(x[!is.na(x)] == TRUE))
data_filtered$losses <- losses

data_filtered$prop_gains <- data_filtered$gains / (data_filtered$gains + data_filtered$losses)

#################################

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_te <- selectSampleType(data_filtered, TE)

#################################

calculatePropMaternal <- function(data, minAge, maxAge) {  
  chroms_affected_frame <- data.frame(matrix(ncol = 4, nrow = 1))
  names(chroms_affected_frame) <- c("maternal_age", "chromosome", "prop", "n")
  age_subset <- data[data$maternal_age > minAge & data$maternal_age <= maxAge & !is.na(data$maternal_age),]
  for (i in 1:23) {
    chrom_subset <- age_subset[age_subset$chroms_affected == i,]
    new_row <- data.frame(cbind(paste(minAge, "-", maxAge), i, mean(chrom_subset$prop_maternal_affected), nrow(chrom_subset)))
    names(new_row) <- c("maternal_age", "chromosome", "prop", "n")
    chroms_affected_frame <- rbind(chroms_affected_frame, new_row)
  }
  chroms_affected_frame <- chroms_affected_frame[-1,]
  chroms_affected_frame$chromosome <- as.numeric(as.character(chroms_affected_frame$chromosome))
  chroms_affected_frame$prop <- as.numeric(as.character(chroms_affected_frame$prop))
  chroms_affected_frame$n <- as.numeric(as.character(chroms_affected_frame$n))
  return(chroms_affected_frame)
}

age18to35 <- calculatePropMaternal(data_blastomere, 18, 35)
age35to50 <- calculatePropMaternal(data_blastomere, 35, 50)
chroms_affected_frame_blastomere <- rbind(age18to35, age35to50)
chroms_affected_frame_blastomere$sample_type <- "Day-3 Blastomere"

age18to35 <- calculatePropMaternal(data_te, 18, 35)
age35to50 <- calculatePropMaternal(data_te, 35, 50)
chroms_affected_frame_te <- rbind(age18to35, age35to50)
chroms_affected_frame_te$sample_type <- "Day-5 TE"

chroms_affected_frame <- rbind(chroms_affected_frame_blastomere, chroms_affected_frame_te)
chroms_affected_frame$se <- sqrt((chroms_affected_frame$prop * (1 - chroms_affected_frame$prop)) / chroms_affected_frame$n)

limits <- aes(ymax = prop + se, ymin = prop - se)

a <- ggplot(data = chroms_affected_frame, aes(x = chromosome, y = prop, color = factor(maternal_age))) + 
  theme_bw() + geom_point() + ylim(0,1) + xlab("Number of Aneuploid Chromosomes") + ylab("\u2640 / (\u2640 + \u2642)") + 
  labs(color = "Maternal Age", linetype="Maternal Age", shape="Maternal Age") + geom_errorbar(limits, width = 0.25) + 
  geom_line() + geom_hline(aes(yintercept = 0.5), lty = "dotted") 

a <- a + facet_grid(. ~ sample_type)

maternal_results <- chroms_affected_frame


#################################################################

calculatePropGains <- function(data, minAge, maxAge) {  
  chroms_affected_frame <- data.frame(matrix(ncol = 4, nrow = 1))
  names(chroms_affected_frame) <- c("maternal_age", "chromosome", "prop", "n")
  age_subset <- data[data$maternal_age > minAge & data$maternal_age <= maxAge & !is.na(data$maternal_age),]
  for (i in 1:23) {
    chrom_subset <- age_subset[age_subset$chroms_affected == i,]
    new_row <- data.frame(cbind(paste(minAge, "-", maxAge), i, mean(chrom_subset$prop_gains), nrow(chrom_subset)))
    names(new_row) <- c("maternal_age", "chromosome", "prop", "n")
    chroms_affected_frame <- rbind(chroms_affected_frame, new_row)
  }
  chroms_affected_frame <- chroms_affected_frame[-1,]
  chroms_affected_frame$chromosome <- as.numeric(as.character(chroms_affected_frame$chromosome))
  chroms_affected_frame$prop <- as.numeric(as.character(chroms_affected_frame$prop))
  chroms_affected_frame$n <- as.numeric(as.character(chroms_affected_frame$n))
  return(chroms_affected_frame)
}

age18to35 <- calculatePropGains(data_blastomere, 18, 35)
age35to50 <- calculatePropGains(data_blastomere, 35, 50)
chroms_affected_frame_blastomere <- rbind(age18to35, age35to50)
chroms_affected_frame_blastomere$sample_type <- "Day-3 Blastomere"

age18to35 <- calculatePropGains(data_te, 18, 35)
age35to50 <- calculatePropGains(data_te, 35, 50)
chroms_affected_frame_te <- rbind(age18to35, age35to50)
chroms_affected_frame_te$sample_type <- "Day-5 TE"

chroms_affected_frame <- rbind(chroms_affected_frame_blastomere, chroms_affected_frame_te)
chroms_affected_frame$se <- sqrt((chroms_affected_frame$prop * (1 - chroms_affected_frame$prop)) / chroms_affected_frame$n)

limits <- aes(ymax = prop + se, ymin = prop - se)

b <- ggplot(data = chroms_affected_frame, aes(x = chromosome, y = prop, color = factor(maternal_age))) + theme_bw() + geom_point() + ylim(0, 1) + xlab("Number of Aneuploid Chromosomes") + ylab("Gains / (Gains + Losses)") + labs(color = "Maternal Age", linetype="Maternal Age", shape="Maternal Age") + geom_errorbar(limits, width = 0.25) + geom_line() + geom_hline(aes(yintercept = 0.5), lty = "dotted")
b <- b + facet_grid(. ~ sample_type)

gains_results <- chroms_affected_frame

####################################

grid.arrange(a, b, nrow = 2)


c <- ggplot(data = chroms_affected_frame, aes(x = chromosome, y = prop, color = factor(maternal_age))) + 
  theme_bw() + ylim(0,1) + xlab("Number of Aneuploid Chromosomes") + ylab("\u2640 / (\u2640 + \u2642)") + 
  labs(color = "Maternal Age", linetype="Maternal Age", shape="Maternal Age")  + 
  geom_blank() + geom_hline(aes(yintercept = 0.5), lty = "dotted") 

c <- c + facet_grid(. ~ sample_type)

d <- ggplot(data = chroms_affected_frame, aes(x = chromosome, y = prop, color = factor(maternal_age))) + geom_blank() + theme_bw() + ylim(0, 1) + xlab("Number of Aneuploid Chromosomes") + ylab("Gains / (Gains + Losses)") + labs(color = "Maternal Age", linetype="Maternal Age", shape="Maternal Age") + geom_hline(aes(yintercept = 0.5), lty = "dotted")
d <- d + facet_grid(. ~ sample_type)

grid.arrange(c, d, nrow = 2)
