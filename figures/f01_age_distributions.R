#################################################################
# This file plots maternal and paternal age distributions.      #
#################################################################

library(ggplot2)
library(gridExtra)
library(gtable)
source("~/Desktop/aneuploidy_analysis-master/aneuploidy_functions.R")

URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" # import the data
url <- getURL(URL)
data <- read.table(textConnection(url), sep=",", header=T)

ages <- data[!duplicated(data$case),][,4:6]

legend_plot <- ggplot(data = ages, aes(x = maternal_age, fill = factor(egg_donor))) + theme_bw() + geom_histogram(position = "dodge", binwidth = 1) + xlab("") + ylab("Count") + scale_x_continuous(limits = c(18, 50)) + scale_fill_discrete(name = "", breaks = c(0, 1), labels = c("Non-donor cases", "Egg donor cases"))

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

legend <- g_legend(legend_plot)

hist_top <- ggplot(data = ages, aes(x = maternal_age, fill = factor(egg_donor))) + theme_bw() + geom_histogram(position = "dodge", binwidth = 1) + xlab("") + ylab("Count") + scale_x_continuous(limits = c(18, 50)) + theme(legend.position="none") 

hist_right <- ggplot(data = ages, aes(x = paternal_age, fill = factor(egg_donor))) + theme_bw() + geom_histogram(position = "dodge", binwidth = 1) + xlab("") + ylab("Count") + scale_x_continuous(limits = c(20, 80)) + coord_flip() + theme(legend.position="none")

scatter <- ggplot(data = ages, aes(x = maternal_age, y = paternal_age, color = factor(egg_donor))) + theme_bw() + geom_point(size = 0.75) + theme(legend.position="none") + xlab('Maternal Age') + ylab('Paternal Age') + scale_y_continuous(limits = c(20, 80)) + scale_x_continuous(limits = c(18, 50))

gA <- ggplot_gtable(ggplot_build(hist_top))
gB <- ggplot_gtable(ggplot_build(hist_right))
gC <- ggplot_gtable(ggplot_build(scatter))
gD <- ggplot_gtable(ggplot_build(legend))

gC$widths <- gA$widths

grid.arrange(gA, legend, gC, gB, nrow = 2, heights = c(2, 4), widths = c(4, 2.2))
grid.arrange(hist_top, legend, scatter, hist_right, ncol=2, nrow=2, widths = c(4, 2.2), heights = c(2, 4))
