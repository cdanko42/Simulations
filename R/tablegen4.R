library(tidyverse)
library(gt)
setwd("..")
auxbias1 <- read.csv("Data/auxcoverage1.csv")[,2:5]
auxbias2 <- read.csv("Data/auxcoverage2.csv")[,2:5]
auxbias3 <- read.csv("Data/auxcoverage3.csv")[,2:5]
auxbias4 <- read.csv("Data/auxcoverage4.csv")[,2:5]
auxbias5 <- read.csv("Data/auxcoverage5.csv")[,2:5]

bias1 <- auxbias1[,1]
bias2 <- auxbias2[, 1]
bias3 <- auxbias3[, 1]
bias4 <- auxbias4[, 1]
bias5 <- auxbias5[, 1]
bias6 <- auxbias1[, 2]
bias7 <- auxbias2[, 2]
bias8 <- auxbias3[, 2]
bias9 <- auxbias4[, 2]
bias10 <- auxbias5[, 2]
bias11 <- auxbias1[, 3]
bias12<- auxbias2[, 3]
bias13 <- auxbias3[, 3]
bias14 <- auxbias4[, 3]
bias15 <- auxbias5[, 3]
bias16 <- auxbias1[, 4]
bias17 <- auxbias2[, 4]
bias18 <- auxbias3[, 4]
bias19 <- auxbias4[, 4]
bias20 <- auxbias5[, 4]


biastable <- as.data.frame(matrix(0, ncol=4, nrow=30))
for (i in 2:6){
  biastable[(6-i)*5+1,1] <- bias1[i]
  biastable[(6-i)*5+2,1] <- bias2[i]
  biastable[(6-i)*5+3,1] <- bias3[i]
  biastable[(6-i)*5+4,1] <- bias4[i]
  biastable[(6-i)*5+5,1] <- bias5[i]
  biastable[(6-i)*5+1,2] <- bias6[i]
  biastable[(6-i)*5+2,2] <- bias7[i]
  biastable[(6-i)*5+3,2] <- bias8[i]
  biastable[(6-i)*5+4,2] <- bias9[i]
  biastable[(6-i)*5+5,2] <- bias10[i]
  biastable[(6-i)*5+1,3] <- bias11[i]
  biastable[(6-i)*5+2,3] <- bias12[i]
  biastable[(6-i)*5+3,3] <- bias13[i]
  biastable[(6-i)*5+4,3] <- bias14[i]
  biastable[(6-i)*5+5,3] <- bias15[i]
  biastable[(6-i)*5+1,4] <- bias16[i]
  biastable[(6-i)*5+2,4] <- bias17[i]
  biastable[(6-i)*5+3,4] <- bias18[i]
  biastable[(6-i)*5+4,4] <- bias19[i]
  biastable[(6-i)*5+5,4] <- bias20[i]
}

colnames(biastable) <- c("2SLPM", "Control Function", "Maximum Likelihood", "Special Regressor")
biastable$names <- 1

biastable$names[1:5] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
biastable$names[6:10] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
biastable$names[11:15] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
biastable$names[16:20] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
biastable$names[21:25] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )

biastable <- biastable[1:25,]
biastable[, 1:4] <- round(biastable[, 1:4], digits=4)
biastable[, 1:4] <- format(biastable[,1:4], scientific=F)
biastable[5, 2] <- NA


gtable <- biastable %>%
  gt(rowname_col= "names") %>%
  tab_header(
    title = "Mean Absolute Deviation for Simulation Estimates"
  ) %>%
tab_row_group(
  group = "Corr(x, µ) =.5",
  rows = 1:5
) %>%
  tab_row_group(
    group = "Corr(x, µ) =.4",
    rows = 6:10
  )%>%
  tab_row_group(
    group = "Corr(x,µ) =.3",
    rows = 11:15
  )%>%
  tab_row_group(
    group = "Corr(x, µ) =.2",
    rows = 16:20
  )%>%
  tab_row_group(
    group = "Corr(x, µ) =.1",
    rows = 21:25
  ) %>%
  tab_options(table.width= pct(66.7))

gtsave(gtable, "Graphics/Table4.pdf")

as_latex(gtable)