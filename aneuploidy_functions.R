#################################################################
# This contains data manipulation functions for working with    #
# aneuploidy call data generated by Natera. 			        #
#################################################################

library(dplyr)
library(data.table)

std <- function(x) sd(x)/sqrt(length(x))

# Function: filterData(data)
# Performs QC filtering to remove samples with excessive low confidence calls
# or whole chromosome nullisomy. Columns that contain aneuploidy calls are 
# hard-coded (7:29, in this case)

filterData <- function(data) {
	data_calls <- data[,7:29]; nrow(data_calls)
	data_calls <- data_calls[(apply(data_calls, 1, function(x) sum(is.na(x))) < 5),] # remove samples with 5+ no-calls
	data_calls <- data_calls[(apply(data_calls, 1, function(x) sum(x[!is.na(x)] == "H000")) != 23),] # remove whole genome nullisomy
	data_filtered <- data[row.names(data_calls),] #recover the original data frame, with all chrom. nullisomy and any chrom. no-calls filtered out
}

filterDataTable <- function(data) {
	data_filtered <- data[apply(data[, 7:29, with = F], 1, function(x) ((sum(is.na(x)) < 5) & (sum(x[!is.na(x)] == "H000") != 23))),]
	return(data_filtered)
}

# Function: selectSampleType(data, sampleType)
# Subsets the dataset based on sample type. 

selectSampleType <- function(data, sampleType) {
	sampleType <- deparse(substitute(sampleType))
	data_subset <- data[data$sample_type == toString(sampleType),]
	return(data_subset)
}

# Function: callPloidy
# Adds a new boolean field called 'ploidy' that indicates whole-chromosome aneuploidy with 'FALSE'
callPloidy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i] != "H110" & data[,i] != "H101" & data[,i+69] != 1)
		aneuploid_frame[,i-6] <- new
	}
	
	chroms_affected <- apply(aneuploid_frame, 1, function(x) sum(x[!is.na(x)]==TRUE))
	data$chroms_affected <- chroms_affected
	
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0) # if any chromosome is aneuploid, call sample aneuploid
	data$ploidy<-TRUE
	data$ploidy[aneuploid_indicator]<-FALSE
	return(data)
}

# Function: callPloidy
# Adds a new boolean field called 'ploidy' that indicates whole-chromosome aneuploidy with 'FALSE'
callPloidy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i] != "H110" & data[,i] != "H101" & data[,i+69] != 1)
		aneuploid_frame[,i-6] <- new
	}
	
	chroms_affected <- apply(aneuploid_frame, 1, function(x) sum(x[!is.na(x)]==TRUE))
	data$chroms_affected <- chroms_affected
	
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0) # if any chromosome is aneuploid, call sample aneuploid
	data$ploidy<-TRUE
	data$ploidy[aneuploid_indicator]<-FALSE
	return(data)
}

# Function: callPloidyTable
# Adds a new boolean field called 'ploidy' that indicates whole-chromosome aneuploidy with 'FALSE'
callPloidyTable <- function(data) {
	data <- data.frame(data)
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i] != "H110" & data[,i] != "H101" & data[,i+69] != 1)
		aneuploid_frame[,i-6] <- new
	}
	
	chroms_affected <- apply(aneuploid_frame, 1, function(x) sum(x[!is.na(x)]==TRUE))
	data$chroms_affected <- chroms_affected
	
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0) # if any chromosome is aneuploid, call sample aneuploid
	data$ploidy<-TRUE
	data$ploidy[aneuploid_indicator]<-FALSE
	return(data.table(data))
}

# Function: callMeiotic
# Adds a new boolean field called 'meiotic' that indicates a BPH aneuploidy (maternal or paternal)
callMeiotic <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i+23] == 1 | data[,i+46] == 1) & (data[,i+69] != 1) # this line for both maternal and paternal
		aneuploid_frame[,i-6] <- new
	}
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0) # if any chromosome is maternalBPH call meiotic
	data$meiotic <- FALSE
	data$meiotic[aneuploid_indicator] <- TRUE
	return(data)
}

# Function: callMaternalMeiotic
# Adds a new boolean field called 'meiotic' that indicates a BPH aneuploidy (maternal)
callMaternalMeiotic <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i+23] == 1) & (data[,i+69] != 1) # this line for maternal
		aneuploid_frame[,i-6] <- new
	}
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0) # if any chromosome is maternalBPH call meiotic
	data$meiotic <- FALSE
	data$meiotic[aneuploid_indicator] <- TRUE
	return(data)
}

# Function: callPaternalMeiotic
# Adds a new boolean field called 'meiotic' that indicates a BPH aneuploidy (paternal)
callPaternalMeiotic <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i+46] == 1) & (data[,i+69] != 1) # this line for paternal
		aneuploid_frame[,i-6] <- new
	}
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0) # if any chromosome is maternalBPH call meiotic
	data$meiotic <- FALSE
	data$meiotic[aneuploid_indicator] <- TRUE
	return(data)
}


# Function: callMitotic
# Adds a new boolean field called 'mitotic' that indicates an aneuploidy affecting a paternal homolog
callMitotic <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i] == "H120" | data[,i] == "H100" | data[,i] == "H102") & (data[,i+46] != 1 & data[,i+69] != 1)
		aneuploid_frame[,i-6] <- new
	}
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0) # if any chromosome is aneuploid, call sample aneuploid
	data$mitotic <- FALSE
	data$mitotic[aneuploid_indicator] <- TRUE
	return(data)
}

# Function: callMaternalTriploidy
# Adds a new boolean field called 'maternalTriploidy' that indicates an aneuploidy affecting a paternal homolog
callMaternalTriploidy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i] == "H210" | data[,i] == "H012") & (data[,i+69] != 1)
		aneuploid_frame[,i-6] <- new
	}
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 19)
	data$maternalTriploidy <- FALSE
	data$maternalTriploidy[aneuploid_indicator] <- TRUE
	return(data)
}


# Function: callEuploidy
# Adds a new boolean field called 'euploidy' that indicates euploidy with 'TRUE'. This is not the 
# inverse of the aneuploid embryos from 'callPloidy', as any segmental duplications or deletions 
# are considered aneuploidies in this case.

callEuploidy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i] != "H110" & data[,i] != "H101")
		aneuploid_frame[,i-6] <- new
	}
	aneuploid_indicator <- (apply(aneuploid_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0) # if any chromosome is aneuploid, call sample aneuploid
	data$euploidy<-TRUE
	data$euploidy[aneuploid_indicator]<-FALSE
	return(data)
}

aneuploidyByAge <- function(data, label) {
	results <- data.frame(matrix(ncol = 5))
	minAge <- min(data$maternal_age[!is.na(data$maternal_age)])
	maxAge <- max(data$maternal_age[!is.na(data$maternal_age)])
	
	for (k in round(minAge):round(maxAge)) {	
		age_subset <- data[round(data$maternal_age) == k & !is.na(round(data$maternal_age)),]
		aneuploidy_prop <- (sum(age_subset$ploidy == FALSE)) / nrow(age_subset)
		se <- sqrt((aneuploidy_prop * (1 - aneuploidy_prop)) / nrow(age_subset)) # calculate standard error of proportion
		aneuploidy <- c(aneuploidy_prop, se, k, label, nrow(age_subset))
		results<-rbind(results, aneuploidy)
	}
	results$X1 <- as.numeric(results$X1)
	results$X2 <- as.numeric(results$X2)
	results$X3 <- as.numeric(results$X3)
	results$X5 <- as.numeric(results$X5)
	return(results[-1,])
}


euploidyByAge <- function(data, label) {
	results <- data.frame(matrix(ncol = 5))
	minAge <- min(data$maternal_age[!is.na(data$maternal_age)])
	maxAge <- max(data$maternal_age[!is.na(data$maternal_age)])
	
	for (k in round(minAge):round(maxAge)) {	
		age_subset <- data[round(data$maternal_age) == k & !is.na(round(data$maternal_age)),]
		euploidy_prop <- (sum(age_subset$euploidy == TRUE)) / nrow(age_subset)
		se <- sqrt((euploidy_prop * (1 - euploidy_prop)) / nrow(age_subset)) # calculate standard error of proportion
		euploidy <- c(euploidy_prop, se, k, label, nrow(age_subset))
		results<-rbind(results, euploidy)
	}
	results$X1 <- as.numeric(results$X1)
	results$X2 <- as.numeric(results$X2)
	results$X3 <- as.numeric(results$X3)
	results$X5 <- as.numeric(results$X5)
	return(results[-1,])
}

meioticByAge <- function(data, label) {
	results <- data.frame(matrix(ncol = 5))
	minAge <- min(data$maternal_age[!is.na(data$maternal_age)])
	maxAge <- max(data$maternal_age[!is.na(data$maternal_age)])
	
	for (k in round(minAge):round(maxAge)) {	
		age_subset <- data[round(data$maternal_age) == k & !is.na(round(data$maternal_age)),]
		aneuploidy_prop <- (sum(age_subset$meiotic == TRUE)) / nrow(age_subset)
		se <- sqrt((aneuploidy_prop * (1 - aneuploidy_prop)) / nrow(age_subset)) # calculate standard error of proportion
		aneuploidy <- c(aneuploidy_prop, se, k, label, nrow(age_subset))
		results<-rbind(results, aneuploidy)
	}
	results$X1 <- as.numeric(results$X1)
	results$X2 <- as.numeric(results$X2)
	results$X3 <- as.numeric(results$X3)
	results$X5 <- as.numeric(results$X5)
	return(results[-1,])
}

paternalMeioticByAge <- function(data, label) {
	results <- data.frame(matrix(ncol = 5))
	minAge <- min(data$paternal_age[!is.na(data$paternal_age)])
	maxAge <- max(data$paternal_age[!is.na(data$paternal_age)])
	
	for (k in round(minAge):round(maxAge)) {	
		age_subset <- data[round(data$paternal_age) == k & !is.na(round(data$paternal_age)),]
		aneuploidy_prop <- (sum(age_subset$meiotic == TRUE)) / nrow(age_subset)
		se <- sqrt((aneuploidy_prop * (1 - aneuploidy_prop)) / nrow(age_subset)) # calculate standard error of proportion
		aneuploidy <- c(aneuploidy_prop, se, k, label, nrow(age_subset))
		results<-rbind(results, aneuploidy)
	}
	results$X1 <- as.numeric(results$X1)
	results$X2 <- as.numeric(results$X2)
	results$X3 <- as.numeric(results$X3)
	results$X5 <- as.numeric(results$X5)
	return(results[-1,])
}


mitoticByAge <- function(data, label) {
	results <- data.frame(matrix(ncol = 5))
	minAge <- min(data$maternal_age[!is.na(data$maternal_age)])
	maxAge <- max(data$maternal_age[!is.na(data$maternal_age)])
	
	for (k in round(minAge):round(maxAge)) {	
		age_subset <- data[round(data$maternal_age) == k & !is.na(round(data$maternal_age)),]
		aneuploidy_prop <- (sum(age_subset$mitotic == TRUE)) / nrow(age_subset)
		se <- sqrt((aneuploidy_prop * (1 - aneuploidy_prop)) / nrow(age_subset)) # calculate standard error of proportion
		aneuploidy <- c(aneuploidy_prop, se, k, label, nrow(age_subset))
		results<-rbind(results, aneuploidy)
	}
	results$X1 <- as.numeric(results$X1)
	results$X2 <- as.numeric(results$X2)
	results$X3 <- as.numeric(results$X3)
	results$X5 <- as.numeric(results$X5)
	return(results[-1,])
}


mitoticByPaternalAge <- function(data, label) {
	results <- data.frame(matrix(ncol = 5))
	minAge <- min(data$paternal_age[!is.na(data$paternal_age)])
	maxAge <- max(data$paternal_age[!is.na(data$paternal_age)])
	
	for (k in round(minAge):round(maxAge)) {	
		age_subset <- data[round(data$paternal_age) == k & !is.na(round(data$paternal_age)),]
		aneuploidy_prop <- (sum(age_subset$mitotic == TRUE)) / nrow(age_subset)
		se <- sqrt((aneuploidy_prop * (1 - aneuploidy_prop)) / nrow(age_subset)) # calculate standard error of proportion
		aneuploidy <- c(aneuploidy_prop, se, k, label, nrow(age_subset))
		results<-rbind(results, aneuploidy)
	}
	results$X1 <- as.numeric(results$X1)
	results$X2 <- as.numeric(results$X2)
	results$X3 <- as.numeric(results$X3)
	results$X5 <- as.numeric(results$X5)
	return(results[-1,])
}


# Function: aneuploidyByCase

aneuploidyByCase <- function(data) {
	data <- data.table(data)
	setorder(data, "case") 
	grouped <- group_by(data, case)
	summary <- summarise(grouped, euploid = sum(ploidy == T), aneuploid = sum(ploidy == F))
	return(summary)
}

# Function: meioticByCase

meioticByCase <- function(data) {
	data <- data.table(data)
	setorder(data, "case") 
	grouped <- group_by(data, case)
	summary <- summarise(grouped, aneuploid = sum(meiotic == T), euploid = sum(meiotic == F))
	return(summary)
}


mitoticByCase <- function(data) {
	data <- data.table(data)
	setorder(data, "case") 
	grouped <- group_by(data, case)
	summary <- summarise(grouped, aneuploid = sum(mitotic == T), euploid = sum(mitotic == F))
	return(summary)
}

# Function: aneuploidyByChromosome
# calculates chromosome-specific proportion of aneuploidy with standard errors 

aneuploidyByChromosome <- function(data) {	
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i] != "H110" & data[,i] != "H101" & data[,i + 69] != 1)
		aneuploid_frame[,i - 6] <- new
	}
	
	p <- apply(aneuploid_frame, 2, function(x) sum(x[!is.na(x)]==TRUE)) / apply(aneuploid_frame, 2, function(x) sum(!is.na(x)))
	n <- apply(aneuploid_frame, 2, function(x) sum(!is.na(x)))

	se <- sqrt((p * (1 - p)) / n)

	results <- data.frame(p, se, c(1:22, "X/Y"))
	names(results) <- c("p", "se", "chrom")
	return(results)
}


# Function: maternalBPHByChromosome
# calculates chromosome-specific proportion of aneuploidy with standard errors 

maternalBPHByChromosome <- function(data) {	
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i+23] == 1) & (data[,i+69] != 1)
		aneuploid_frame[,i - 6] <- new
	}
	
	p <- apply(aneuploid_frame, 2, function(x) sum(x[!is.na(x)]==TRUE)) / apply(aneuploid_frame, 2, function(x) sum(!is.na(x)))
	n <- apply(aneuploid_frame, 2, function(x) sum(!is.na(x)))

	se <- sqrt((p * (1 - p)) / n)

	results <- data.frame(p, se, c(1:22, "X"))
	names(results) <- c("p", "se", "chrom")
	return(results)
}

# Function: paternalBPHByChromosome
# calculates chromosome-specific proportion of aneuploidy with standard errors 

paternalBPHByChromosome <- function(data) {	
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i+46] == 1) & (data[,i+69] != 1)
		aneuploid_frame[,i - 6] <- new
	}
	
	p <- apply(aneuploid_frame, 2, function(x) sum(x[!is.na(x)]==TRUE)) / apply(aneuploid_frame, 2, function(x) sum(!is.na(x)))
	n <- apply(aneuploid_frame, 2, function(x) sum(!is.na(x)))

	se <- sqrt((p * (1 - p)) / n)

	results <- data.frame(p, se, c(1:22, "Y"))
	names(results) <- c("p", "se", "chrom")
	return(results)
}


# Function: mitoticByChromosome
# calculates chromosome-specific proportion of aneuploidy with standard errors 

mitoticByChromosome <- function(data) {	
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[,i] == "H120" | data[,i] == "H100" | data[,i] == "H102") & (data[,i+46] != 1 & data[,i+69] != 1)
		aneuploid_frame[,i - 6] <- new
	}
	
	p <- apply(aneuploid_frame, 2, function(x) sum(x[!is.na(x)]==TRUE)) / apply(aneuploid_frame, 2, function(x) sum(!is.na(x)))
	n <- apply(aneuploid_frame, 2, function(x) sum(!is.na(x)))

	se <- sqrt((p * (1 - p)) / n)

	results <- data.frame(p, se, c(1:22, "X/Y"))
	names(results) <- c("p", "se", "chrom")
	return(results)
}
