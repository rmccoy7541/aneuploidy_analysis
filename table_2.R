library(ggplot2)
library(gridExtra)
library(RCurl)
library(data.table)

source('~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R', chdir = TRUE)

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- fread(url, sep=",", header=T)

data_filtered <- filterDataTable(data)

data_filtered <- callPloidyTable(data_filtered)

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

maternalErrs <- function(data) {
	maternal_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
	  new <- (data[, i, with = F] == "H200" | data[, i, with = F] == "H020" | data[, i, with = F] == "H010" | data[, i, with = F] == "H001" | data[, i, with = F] == "H000" | data[, i, with = F] == "H210" | data[, i, with = F] == "H201" | data[, i, with = F] == "H021") & (data[, i + 69, with = F] != 1)
	  maternal_frame[, i - 6] <- new
	}
	return(rowSums(maternal_frame, na.rm = T))
}

paternalErrs <- function(data) {
	paternal_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
	  new <- (data[, i, with = F] == "H200" | data[, i, with = F] == "H020" | data[, i, with = F] == "H100" | data[, i, with = F] == "H000" | data[, i, with = F] == "H102" | data[, i, with = F] == "H120" | data[, i, with = F] == "H201" | data[, i, with = F] == "H021" | data[, i, with = F] == "H111") & (data[, i + 69, with = F] != 1)
	  paternal_frame[,i - 6] <- new
	}
	return(rowSums(paternal_frame, na.rm = T))
}


totalChroms <- function(data) {
	chroms_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	aneuploid_frame <- data[, 7:29, with = F]
	chroms_frame[aneuploid_frame == "H110"] <- 2
	chroms_frame[aneuploid_frame == "H101"] <- 2
	chroms_frame[aneuploid_frame == "H011"] <- 2
	chroms_frame[aneuploid_frame == "H210"] <- 3
	chroms_frame[aneuploid_frame == "H120"] <- 3
	chroms_frame[aneuploid_frame == "H111"] <- 3
	chroms_frame[aneuploid_frame == "H201"] <- 3
	chroms_frame[aneuploid_frame == "H102"] <- 3
	chroms_frame[aneuploid_frame == "H100"] <- 1
	chroms_frame[aneuploid_frame == "H010"] <- 1
	chroms_frame[aneuploid_frame == "H001"] <- 1
	chroms_frame[aneuploid_frame == "H000"] <- 0
	return(rowSums(chroms_frame, na.rm = T))
}

totalMatChroms <- function(data) {
	chroms_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	aneuploid_frame <- data[, 7:29, with = F]
	chroms_frame[aneuploid_frame == "H110"] <- 1
	chroms_frame[aneuploid_frame == "H101"] <- 1
	chroms_frame[aneuploid_frame == "H011"] <- 0
	chroms_frame[aneuploid_frame == "H210"] <- 2
	chroms_frame[aneuploid_frame == "H120"] <- 1
	chroms_frame[aneuploid_frame == "H111"] <- 1
	chroms_frame[aneuploid_frame == "H201"] <- 2
	chroms_frame[aneuploid_frame == "H102"] <- 1
	chroms_frame[aneuploid_frame == "H100"] <- 1
	chroms_frame[aneuploid_frame == "H010"] <- 0
	chroms_frame[aneuploid_frame == "H001"] <- 0
	chroms_frame[aneuploid_frame == "H000"] <- 0
	return(rowSums(chroms_frame, na.rm = T))
}

totalPatChroms <- function(data) {
	chroms_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	aneuploid_frame <- data[, 7:29, with = F]
	chroms_frame[aneuploid_frame == "H110"] <- 1
	chroms_frame[aneuploid_frame == "H101"] <- 1
	chroms_frame[aneuploid_frame == "H011"] <- 2
	chroms_frame[aneuploid_frame == "H210"] <- 1
	chroms_frame[aneuploid_frame == "H120"] <- 2
	chroms_frame[aneuploid_frame == "H111"] <- 2
	chroms_frame[aneuploid_frame == "H201"] <- 1
	chroms_frame[aneuploid_frame == "H102"] <- 2
	chroms_frame[aneuploid_frame == "H100"] <- 0
	chroms_frame[aneuploid_frame == "H010"] <- 1
	chroms_frame[aneuploid_frame == "H001"] <- 1
	chroms_frame[aneuploid_frame == "H000"] <- 0
	return(rowSums(chroms_frame, na.rm = T))
}

data_blastomere$maternalChroms <- totalMatChroms(data_blastomere)
data_blastomere$paternalChroms <- totalPatChroms(data_blastomere)
data_blastomere$totalChroms <- totalChroms(data_blastomere)
data_te$maternalChroms <- totalMatChroms(data_te)
data_te$paternalChroms <- totalPatChroms(data_te)
data_te$totalChroms <- totalChroms(data_te)

set.seed(42)
data_sampled <- rbind(data_blastomere[sample(nrow(data_te)),], data_te)

df <- data.frame(mat = data_sampled$maternalChroms, pat = data_sampled$paternalChroms, sample_type = data_sampled$sample_type)
levels(df$sample_type) <- c("Day-3 Blastomere", "Day-5 TE Biopsy")
p <- ggplot(df, aes(x = mat, y = pat)) + stat_binhex() + scale_fill_gradientn(colours = rev(rainbow(3)), name = "Samples", trans = "log", breaks = 10^(0:6))
p + facet_grid(. ~ sample_type) + theme_bw() + ylab('Number of Paternal Chromosomes') + xlab('Number of Maternal Chromosomes')

####################################################

trisomy <- function(data) {
	aneuploid_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data)))
	for (i in 7:29) {
		new <- ((data[, i, with = F] == "H120") | (data[, i, with = F] == "H102") | (data[, i, with = F] == "H210")) & (data[, i + 69, with = F] != 1)
		aneuploid_frame[, i - 6] <- new
	}
	return(rowSums(aneuploid_frame, na.rm = T))
}

sum(trisomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) 
singleTrisomy <- sum(trisomy(data_blastomere) == 1 & aneuploidChroms(data_blastomere) == 1) / nrow(data_blastomere)
se(singleTrisomy, nrow(data_blastomere))

sum(trisomy(data_te) == 1 & aneuploidChroms(data_te) == 1) 
singleTrisomy <- sum(trisomy(data_te) == 1 & aneuploidChroms(data_te) == 1) / nrow(data_te)
se(singleTrisomy, nrow(data_te))

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

sum(trisomy(data_blastomere) > 19)
triploidy <- sum(trisomy(data_blastomere) > 19) / nrow(data_blastomere)
se(triploidy, nrow(data_blastomere))

sum(trisomy(data_te) > 19)
triploidy <- sum(trisomy(data_te) > 19) / nrow(data_te)
se(triploidy, nrow(data_te))

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

sum(monosomy(data_blastomere) > 19)
haploidy <- sum(monosomy(data_blastomere) > 19) / nrow(data_blastomere)
se(haploidy, nrow(data_blastomere))

sum(monosomy(data_te) > 19)
haploidy <- sum(monosomy(data_te) > 19) / nrow(data_te)
se(haploidy, nrow(data_te))

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

####################################################

sum(aneuploidChroms(data_blastomere) > 2 & aneuploidChroms(data_blastomere) < 20 )
complex <- sum(aneuploidChroms(data_blastomere) > 2 & aneuploidChroms(data_blastomere) < 20 ) / nrow(data_blastomere)
se(complex, nrow(data_blastomere))

sum(aneuploidChroms(data_te) > 2 & aneuploidChroms(data_te) < 20 )
complex <- sum(aneuploidChroms(data_te) > 2 & aneuploidChroms(data_te) < 20 ) / nrow(data_te)
se(complex, nrow(data_te))
