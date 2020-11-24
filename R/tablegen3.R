library(tidyverse)
library(gt)
setwd("..")
endogeneity1 <- read.csv("Data/endo1.csv")[,2]
endogeneity2 <- read.csv("Data/endo2.csv")[,2]
endogeneity3 <- read.csv("Data/endo3.csv")[,2]
endogeneity4 <- read.csv("Data/endo4.csv")[,2]
endogeneity5 <- read.csv("Data/endo5.csv")[,2]
endogeneity6 <- read.csv("Data/endo6.csv")[,2]
endogeneity7 <- read.csv("Data/endo7.csv")[,2]
endogeneity8 <- read.csv("Data/endo8.csv")[,2]
endogeneity9 <- read.csv("Data/endo9.csv")[,2]
endogeneity10 <- read.csv("Data/endo10.csv")[,2]
endogeneity11 <- read.csv("Data/endo11.csv")[,2]
endogeneity12 <- read.csv("Data/endo12.csv")[,2]
endogeneity13 <- read.csv("Data/endo13.csv")[,2]
endogeneity14 <- read.csv("Data/endo14.csv")[,2]
endogeneity15 <- read.csv("Data/endo15.csv")[,2]
endogeneity16 <- read.csv("Data/endogeneity16.csv")[,2]
endogeneity17 <- read.csv("Data/endogeneity17.csv")[,2]
endogeneity18 <- read.csv("Data/endogeneity18.csv")[,2]
endogeneity19 <- read.csv("Data/endogeneity19.csv")[,2]
endogeneity20 <- read.csv("Data/endogeneity20.csv")[,2]


endogeneitytable <- as.data.frame(matrix(0, ncol=4, nrow=30))
for (i in 2:6){
  endogeneitytable[(6-i)*5+1,1] <- endogeneity1[i]/1000
  endogeneitytable[(6-i)*5+2,1] <- endogeneity2[i]/1000
  endogeneitytable[(6-i)*5+3,1] <- endogeneity3[i]/1000
  endogeneitytable[(6-i)*5+4,1] <- endogeneity4[i]/1000
  endogeneitytable[(6-i)*5+5,1] <- endogeneity5[i]/1000
  endogeneitytable[(6-i)*5+1,2] <- endogeneity6[i]/1000
  endogeneitytable[(6-i)*5+2,2] <- endogeneity7[i]/1000
  endogeneitytable[(6-i)*5+3,2] <- endogeneity8[i]/1000
  endogeneitytable[(6-i)*5+4,2] <- endogeneity9[i]/1000
  endogeneitytable[(6-i)*5+5,2] <- endogeneity10[i]/1000
  endogeneitytable[(6-i)*5+1,3] <- endogeneity11[i]/1000
  endogeneitytable[(6-i)*5+2,3] <- endogeneity12[i]/1000
  endogeneitytable[(6-i)*5+3,3] <- endogeneity13[i]/1000
  endogeneitytable[(6-i)*5+4,3] <- endogeneity14[i]/1000
  endogeneitytable[(6-i)*5+5,3] <- endogeneity15[i]/1000
  endogeneitytable[(6-i)*5+1,4] <- endogeneity16[i]/1000
  endogeneitytable[(6-i)*5+2,4] <- endogeneity17[i]/1000
  endogeneitytable[(6-i)*5+3,4] <- endogeneity18[i]/1000
  endogeneitytable[(6-i)*5+4,4] <- endogeneity19[i]/1000
  endogeneitytable[(6-i)*5+5,4] <- endogeneity20[i]/1000
}

colnames(endogeneitytable) <- c("2SLPM", "Control Function", "Maximum Likelihood", "Special Regressor")
endogeneitytable$names <- 1

endogeneitytable$names[1:5] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
endogeneitytable$names[6:10] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
endogeneitytable$names[11:15] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
endogeneitytable$names[16:20] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )
endogeneitytable$names[21:25] <- c("Corr(x, z) = .1", "Corr(x, z) = .3", "Corr(x, z) = .5", "Corr(x, z) = .7", "Corr(x, z) = .9"  )

endogeneitytable <- endogeneitytable[1:25,]
endogeneitytable[, 1:4] <- round(endogeneitytable[, 1:4], digits=4)
endogeneitytable[, 1:4] <- format(endogeneitytable[,1:4], scientific=F)
endogeneitytable[5, 2] <- NA


gtable <- endogeneitytable %>%
  gt(rowname_col= "names") %>%
  tab_header(
    title = "Results from Simulated Endogeneity Tests"
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

gtsave(gtable, "Graphics/Table3.pdf")

as_latex(gtable)