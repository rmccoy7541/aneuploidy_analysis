#################################################################
# This file generates correlation plots of chromosome length    #
# and chromosome-specific rates of meiotic- and mitotic-origin  #
# aneuploidy.                                                   #
#################################################################

library(ggplot2)
library(gridExtra)
library(gtable)
library(RCurl)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

data_filtered <- filterData(data)

#################################################################

data_blastomere <- selectSampleType(data_filtered, blastomere)
data_te <- selectSampleType(data_filtered, TE)

blastomere_meiotic <- data.frame(matrix(ncol = 23, nrow = nrow(data_blastomere)))
for (i in 7:29) {
	new <- (data_blastomere[,i+23] == 1) & (data_blastomere[,i+69] != 1) # this line for maternal
  	blastomere_meiotic[,i - 6] <- new
}

te_meiotic <- data.frame(matrix(ncol = 23, nrow = nrow(data_te)))
for (i in 7:29) {
	new <- (data_te[,i+23] == 1) & (data_te[,i+69] != 1) # this line for maternal
  	te_meiotic[,i - 6] <- new
}

blastomere_mitotic <- data.frame(matrix(ncol = 23, nrow = nrow(data_blastomere)))
for (i in 7:29) {
	new <- (data_blastomere[,i] == "H120" | data_blastomere[,i] == "H100" | data_blastomere[,i] == "H102") & (data_blastomere[,i+46] != 1 & data_blastomere[,i+69] != 1)
  	blastomere_mitotic[,i - 6] <- new
}

te_mitotic <- data.frame(matrix(ncol = 23, nrow = nrow(data_te)))
for (i in 7:29) {
	new <- (data_te[,i] == "H120" | data_te[,i] == "H100" | data_te[,i] == "H102") & (data_te[,i+46] != 1 & data_te[,i+69] != 1)
  	te_mitotic[,i - 6] <- new
}

lengthCorrelation <- function(data_frame, sexChrom) {
  chromLengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566)
  n <- apply(data_frame, 2, function(x) sum(!is.na(x)))
  if (sexChrom == "X") {
    chromLengths <- c(chromLengths, 155270560)
  } else if (sexChrom == "Y") {
    chromLengths <- c(chromLengths, 59373566)
  }
  lengths <- data.frame(apply(data_frame, 2, function(x) sum(x[!is.na(x)] == TRUE)),apply(data_frame, 2, function(x) sum(!is.na(x))), chromLengths, c(1:22, sexChrom))
  names(lengths) <- c("affected", "total","length", "chrom")
  
  print(cor.test(lengths$affected / lengths$total, lengths$length))
  
  prop <- lengths$affected / lengths$total
  se <- sqrt((prop * (1 - prop)) / n)
  lengths$prop <- prop
  lengths$se <- se
  return(lengths)
}

meiotic_blastomere_lengths <- lengthCorrelation(blastomere_meiotic, "X") 
meiotic_blastomere_lengths$type <- "Day-3 Blastomere"
meiotic_te_lengths <- lengthCorrelation(te_meiotic, "X") 
meiotic_te_lengths$type <- "Day-5 TE"
meiotic_lengths <- rbind(meiotic_blastomere_lengths, meiotic_te_lengths)

limits <- aes(ymax = prop + se, ymin = prop - se, xmin = length, xmax = length)

a <- ggplot(data = meiotic_lengths, aes(y = prop, x = length)) + theme_bw() + geom_errorbar(limits) + theme(axis.text.x = element_text(angle = 65, hjust = 1)) + ylab("Prop. with Maternal BPH Error") + xlab("Chromosome Length (bp)") + geom_text(aes(label = chrom), color = "gray")

a + facet_grid(. ~ type)

mitotic_blastomere_lengths <- lengthCorrelation(blastomere_mitotic, "Y") 
mitotic_blastomere_lengths$type <- "Day-3 Blastomere"
mitotic_te_lengths <- lengthCorrelation(te_mitotic, "Y") 
mitotic_te_lengths$type <- "Day-5 TE"
mitotic_lengths <- rbind(mitotic_blastomere_lengths, mitotic_te_lengths)

limits <- aes(ymax = prop + se, ymin = prop - se, xmin = length, xmax = length)

b <- ggplot(data = mitotic_lengths, aes(y = prop, x = length)) + theme_bw() + geom_errorbar(limits) + theme(axis.text.x = element_text(angle = 65, hjust = 1)) + ylab("Prop. with Mitotic Error") + xlab("Chromosome Length (bp)") + geom_text(aes(label = chrom), color = "gray")

b + facet_grid(. ~ type)

cor.test(meiotic_blastomere_lengths[c(1:14, 17:20, 23),]$prop, meiotic_blastomere_lengths[c(1:14, 17:20, 23),]$length) 

cor.test(meiotic_te_lengths[c(1:14, 17:20, 23),]$prop, meiotic_te_lengths[c(1:14, 17:20, 23),]$length) 
