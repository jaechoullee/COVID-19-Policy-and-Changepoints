
library(tidyverse)
library(readxl)
library(fitdistrplus)

########################################################################################################
# CASES #####
########################################################################################################


##############################################################
##########################   X1  ####################################
##############################################################

data.X1.Cases <-read_xlsx("Random_Variable_Confirmed_Cases.xlsx", sheet = "X1")

data.X1.Cases <- c(data.X1.Cases$diff1,data.X1.Cases$diff2,data.X1.Cases$diff3,data.X1.Cases$diff4,data.X1.Cases$diff5,
                   data.X1.Cases$diff6,data.X1.Cases$diff7,data.X1.Cases$diff8,data.X1.Cases$diff9,data.X1.Cases$diff10,
                   data.X1.Cases$diff11,data.X1.Cases$diff12)

data.X1.Cases<-data.X1.Cases[!is.na(data.X1.Cases)]

mean(data.X1.Cases,na.rm = TRUE)


# hypothesis test do 3, 7, 10, 14

wilcox.test(data.X1.Cases, mu = 3, alternative = "greater") 

wilcox.test(data.X1.Cases, mu = 7, alternative = "greater") 

wilcox.test(data.X1.Cases, mu = 10, alternative = "greater") 

wilcox.test(data.X1.Cases, mu = 14, alternative = "greater") 


fitdist.X1.cases <- fitdist(data.X1.Cases, "lnorm")
summary(fitdist.X1.cases)

par(mfrow = c(2, 2))

denscomp(list(fitdist.X1.cases))
qqcomp(list(fitdist.X1.cases))
cdfcomp(list(fitdist.X1.cases))
ppcomp(list(fitdist.X1.cases))

gofstat(fitdist.X1.cases)

fitdist.X1.cases.weibull <- fitdist(data.X1.Cases, "weibull")
fitdist.X1.cases.gamma <- fitdist(data.X1.Cases, "gamma")


dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
fitdist.X1.cases.gumbel <- fitdist(data.X1.Cases,  start=list(a=10, b=10), "gumbel")
fitdist.X1.cases.lognormal <- fitdist(data.X1.Cases, "lnorm")

gofstat(list(fitdist.X1.cases.weibull, fitdist.X1.cases.gamma, fitdist.X1.cases.gumbel,
             fitdist.X1.cases.lognormal), fitnames = c("Weibull", "Gamma", "Gumbel", "Lognormal"))




########################################################################################################
# DEATHS #####
########################################################################################################

##############################################################
##########################   X1  ####################################
##############################################################

data.X1.Deaths <-read_xlsx("Random_Variable_Deaths.xlsx", sheet = "X1")

data.X1.Deaths <- c(data.X1.Deaths$diff1,data.X1.Deaths$diff2,data.X1.Deaths$diff3,data.X1.Deaths$diff4,data.X1.Deaths$diff5,data.X1.Deaths$diff6,data.X1.Deaths$diff7,
                    data.X1.Deaths$diff8,data.X1.Deaths$diff9,data.X1.Deaths$diff10,data.X1.Deaths$diff11)

data.X1.Deaths<-data.X1.Deaths[!is.na(data.X1.Deaths)]

mean(data.X1.Deaths,na.rm = TRUE)


# hypothesis test do 3, 7, 10, 14

wilcox.test(data.X1.Deaths, mu = 3, alternative = "greater") 

wilcox.test(data.X1.Deaths, mu = 7, alternative = "greater") 

wilcox.test(data.X1.Deaths, mu = 10, alternative = "greater") 

wilcox.test(data.X1.Deaths, mu = 14, alternative = "greater") 

fitdist.X1.Deaths <- fitdist(data.X1.Deaths, "lnorm")

summary(fitdist.X1.Deaths)

par(mfrow = c(2, 2))

denscomp(list(fitdist.X1.Deaths))
qqcomp(list(fitdist.X1.Deaths))
cdfcomp(list(fitdist.X1.Deaths))
ppcomp(list(fitdist.X1.Deaths))

gofstat(fitdist.X1.Deaths)

fitdist.X1.deaths.weibull <- fitdist(data.X1.Deaths, "weibull")
fitdist.X1.deaths.gamma <- fitdist(data.X1.Deaths, "gamma")

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
fitdist.X1.deaths.gumbel <- fitdist(data.X1.Deaths,  start=list(a=10, b=10), "gumbel")
fitdist.X1.deaths.lognormal <- fitdist(data.X1.Deaths, "lnorm")

gofstat(list(fitdist.X1.deaths.weibull, fitdist.X1.deaths.gamma, fitdist.X1.deaths.gumbel, fitdist.X1.deaths.lognormal), fitnames = c("Weibull", "Gamma", "Gumbel", "Lognormal"))



