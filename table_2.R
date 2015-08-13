library(ggplot2)
library(gridExtra)
library(RCurl)
library(data.table)

source('~/Desktop/aneuploidy-analysis/aneuploidy_functions.R', chdir = TRUE)

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- fread(url, sep=",", header=T)

data_filtered <- filterDataTable(data)

data_blastomere <- selectSampleType(data_filtered, blastomere)
data_te <- selectSampleType(data_filtered, TE)

####################################################

se <- function(p, n) {
	sqrt((p * (1 - p)) / n)
}

####################################################

aneuploidChroms <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[, i, with = F] != "H110") & (data[, i, with = F] != "H101") & (data[, i + 69, with = F] != 1) & !is.na(data[, i, with = F])
		aneuploid_frame[, i - 6] <- new
	}
	return(rowSums(aneuploid_frame, na.rm = T))
}

####################################################

maternalTrisomy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[, i, with = F] == "H210") & (data[, i + 69, with = F] != 1)
		aneuploid_frame[, i - 6] <- new
	}
	return(rowSums(aneuploid_frame, na.rm = T))
}

sum(maternalTrisomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) 
singleMatTrisomy <- sum(maternalTrisomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) / nrow(data_blastomere)
se(singleMatTrisomy, nrow(data_blastomere))

sum(maternalTrisomy(data_te) == 1 & aneuploidChroms(data_te) == 1) 
singleMatTrisomy <- sum(maternalTrisomy(data_te) == 1 & aneuploidChroms(data_te) == 1) / nrow(data_te)
se(singleMatTrisomy, nrow(data_te))

####################################################

paternalTrisomy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- ((data[, i, with = F] == "H120") | (data[, i, with = F] == "H102")) & (data[, i + 69, with = F] != 1)
		aneuploid_frame[, i - 6] <- new
	}
	return(rowSums(aneuploid_frame, na.rm = T))
}

sum(paternalTrisomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) 
singlePatTrisomy <- sum(paternalTrisomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) / nrow(data_blastomere)
se(singlePatTrisomy, nrow(data_blastomere))

sum(paternalTrisomy(data_te) == 1 & aneuploidChroms(data_te) == 1) 
singlePatTrisomy <- sum(paternalTrisomy(data_te) == 1 & aneuploidChroms(data_te) == 1) / nrow(data_te)
se(singlePatTrisomy, nrow(data_te))

####################################################

maternalMonosomy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[, i, with = F] == "H010") & (data[, i + 69, with = F] != 1)
		aneuploid_frame[, i-6] <- new
	}
	return(rowSums(aneuploid_frame, na.rm = T))
}

sum(maternalMonosomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) 
singleMatMonosomy <- sum(maternalMonosomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) / nrow(data_blastomere)
se(singleMatMonosomy, nrow(data_blastomere))

sum(maternalMonosomy(data_te) == 1 & aneuploidChroms(data_te) == 1) 
singleMatMonosomy <- sum(maternalMonosomy(data_te) == 1 & aneuploidChroms(data_te) == 1) / nrow(data_te)
se(singleMatMonosomy, nrow(data_te))

####################################################

paternalMonosomy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[, i, with = F] == "H100") & (data[, i + 69, with = F] != 1)
		aneuploid_frame[, i - 6] <- new
	}
	return(rowSums(aneuploid_frame, na.rm = T))
}

sum(paternalMonosomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) 
singlePatMonosomy <- sum(paternalMonosomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) / nrow(data_blastomere)
se(singlePatMonosomy, nrow(data_blastomere))

sum(paternalMonosomy(data_te) == 1 & aneuploidChroms(data_te) == 1) 
singlePatMonosomy <- sum(paternalMonosomy(data_te) == 1 & aneuploidChroms(data_te) == 1) / nrow(data_te)
se(singlePatMonosomy, nrow(data_te))

####################################################

nullisomy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- (data[, i, with = F] == "H000") & (data[, i + 69, with = F] != 1)
		aneuploid_frame[, i - 6] <- new
	}
	return(rowSums(aneuploid_frame, na.rm = T))
}

sum(nullisomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) 
singleNullisomy <- sum(nullisomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) / nrow(data_blastomere)
se(singleNullisomy, nrow(data_blastomere))

sum(nullisomy(data_te) == 1 & aneuploidChroms(data_te) == 1) 
singleNullisomy <- sum(nullisomy(data_te) == 1 & aneuploidChroms(data_te) == 1) / nrow(data_te)
se(singleNullisomy, nrow(data_te))

####################################################

sum(maternalTrisomy(data_blastomere) > 19)
matTriploidy <- sum(maternalTrisomy(data_blastomere) > 19) / nrow(data_blastomere)
se(matTriploidy, nrow(data_blastomere))

sum(maternalTrisomy(data_te) > 19)
matTriploidy <- sum(maternalTrisomy(data_te) > 19) / nrow(data_te)
se(matTriploidy, nrow(data_te))

####################################################

sum(paternalTrisomy(data_blastomere) > 19)
patTriploidy <- sum(paternalTrisomy(data_blastomere) > 19) / nrow(data_blastomere)
se(patTriploidy, nrow(data_blastomere))

sum(paternalTrisomy(data_te) > 19)
patTriploidy <- sum(paternalTrisomy(data_te) > 19) / nrow(data_te)
se(patTriploidy, nrow(data_te))

####################################################

sum(maternalMonosomy(data_blastomere) > 19)
matHaploidy <- sum(maternalMonosomy(data_blastomere) > 19) / nrow(data_blastomere)
se(matHaploidy, nrow(data_blastomere))

sum(maternalMonosomy(data_te) > 19)
matHaploidy <- sum(maternalMonosomy(data_te) > 19) / nrow(data_te)
se(matHaploidy, nrow(data_te))

####################################################

sum(paternalMonosomy(data_blastomere) > 19)
patHaploidy <- sum(paternalMonosomy(data_blastomere) > 19) / nrow(data_blastomere)
se(patHaploidy, nrow(data_blastomere))

sum(paternalMonosomy(data_te) > 19)
patHaploidy <- sum(paternalMonosomy(data_te) > 19) / nrow(data_te)
se(patHaploidy, nrow(data_te))

