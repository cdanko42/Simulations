library(tidyverse)
library(gt)
setwd("..")
coverage1 <- read.csv("Data/coverage1.csv")[,2]
coverage2 <- read.csv("Data/coverage2.csv")[,2]
coverage3 <- read.csv("Data/coverage3.csv")[,2]
coverage4 <- read.csv("Data/coverage4.csv")[,2]
coverage5 <- read.csv("Data/coverage5.csv")[,2]
coverage6 <- read.csv("Data/coverage6.csv")[,2]
coverage7 <- read.csv("Data/coverage7.csv")[,2]
coverage8 <- read.csv("Data/coverage8.csv")[,2]
coverage9 <- read.csv("Data/coverage9.csv")[,2]
coverage10 <- read.csv("Data/coverage10.csv")[,2]
coverage11 <- read.csv("Data/coverage11.csv")[,2]
coverage12 <- read.csv("Data/coverage12.csv")[,2]
coverage13 <- read.csv("Data/coverage13.csv")[,2]
coverage14 <- read.csv("Data/coverage14.csv")[,2]
coverage15 <- read.csv("Data/coverage15.csv")[,2]
coverage16 <- read.csv("Data/coverage16.csv")[,2]
coverage17 <- read.csv("Data/coverage17.csv")[,2]
coverage18 <- read.csv("Data/coverage18.csv")[,2]
coverage19 <- read.csv("Data/coverage18.csv")[,2]
coverage20 <- read.csv("Data/coverage20.csv")[,2]


coveragetable <- as.data.frame(matrix(0, ncol=4, nrow=30))
for (i in 2:6){
  coveragetable[(6-i)*5+1,1] <- coverage1[i]/1000
  coveragetable[(6-i)*5+2,1] <- coverage2[i]/1000
  coveragetable[(6-i)*5+3,1] <- coverage3[i]/1000
  coveragetable[(6-i)*5+4,1] <- coverage4[i]/1000
  coveragetable[(6-i)*5+5,1] <- coverage5[i]/1000
  coveragetable[(6-i)*5+1,2] <- coverage6[i]/1000
  coveragetable[(6-i)*5+2,2] <- coverage7[i]/1000
  coveragetable[(6-i)*5+3,2] <- coverage8[i]/1000
  coveragetable[(6-i)*5+4,2] <- coverage9[i]/1000
  coveragetable[(6-i)*5+5,2] <- coverage10[i]/1000
  coveragetable[(6-i)*5+1,3] <- coverage11[i]/1000
  coveragetable[(6-i)*5+2,3] <- coverage12[i]/1000
  coveragetable[(6-i)*5+3,3] <- coverage13[i]/1000
  coveragetable[(6-i)*5+4,3] <- coverage14[i]/1000
  coveragetable[(6-i)*5+5,3] <- coverage15[i]/1000
  coveragetable[(6-i)*5+1,4] <- coverage16[i]/1000
  coveragetable[(6-i)*5+2,4] <- coverage17[i]/1000
  coveragetable[(6-i)*5+3,4] <- coverage18[i]/1000
  coveragetable[(6-i)*5+4,4] <- coverage19[i]/1000
  coveragetable[(6-i)*5+5,4] <- coverage20[i]/1000
}

colnames(coveragetable) <- c("2SLPM", "Control Function", "Maximum Likelihood", "Special Regressor")
coveragetable$names <- 1

coveragetable$names[1:5] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
coveragetable$names[6:10] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
coveragetable$names[11:15] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
coveragetable$names[16:20] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
coveragetable$names[21:25] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )

coveragetable <- coveragetable[1:25,]
coveragetable[, 1:4] <- round(coveragetable[, 1:4], digits=4)
coveragetable[, 1:4] <- format(coveragetable[,1:4], scientific=F)
coveragetable[5, 2] <- NA


gtable <- coveragetable %>%
  gt(rowname_col= "names") %>%
  tab_header(
    title = "Coverage of Simulation Estimates"
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

gtsave(gtable, "Graphics/Table2.pdf")

as_latex(gtable)