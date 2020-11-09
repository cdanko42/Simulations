library(tidyverse)
library(gt)
setwd("..")
bias1 <- read.csv("Data/bias1.csv")[,2]
bias2 <- read.csv("Data/bias2.csv")[,2]
bias3 <- read.csv("Data/bias3.csv")[,2]
bias4 <- read.csv("Data/bias4.csv")[,2]
bias5 <- read.csv("Data/bias5.csv")[,2]
bias6 <- read.csv("Data/bias6.csv")[,2]
bias7 <- read.csv("Data/bias7.csv")[,2]
bias8 <- read.csv("Data/bias8.csv")[,2]
bias9 <- read.csv("Data/bias9.csv")[,2]
bias10 <- read.csv("Data/bias10.csv")[,2]
bias11 <- read.csv("Data/bias11.csv")[,2]
bias12 <- read.csv("Data/bias12.csv")[,2]
bias13 <- read.csv("Data/bias13.csv")[,2]
bias14 <- read.csv("Data/bias14.csv")[,2]
bias15 <- read.csv("Data/bias15.csv")[,2]
bias16 <- read.csv("Data/bias16.csv")[,2]
bias17 <- read.csv("Data/bias17.csv")[,2]
bias18 <- read.csv("Data/bias18.csv")[,2]
bias19 <- c(0,0,0,0,0,0)
bias20 <- read.csv("Data/bias20.csv")[,2]


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

colnames(biastable) <- c("2SLPM", "Control Function", "MLE", "Special Regressor")
biastable$names <- 1

biastable$names[1:5] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
biastable$names[6:10] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
biastable$names[11:15] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
biastable$names[16:20] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
biastable$names[21:25] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )

biastable <- biastable[1:25,]
biastable[, 1:4] <- round(biastable[, 1:4], digits=4)

gtable <- biastable %>%
  gt(rowname_col= "names") %>%
  tab_header(
    title = "Summary of Binomial Tests"
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
  )

gtsave(gtable, "Graphics/Table1.png")

