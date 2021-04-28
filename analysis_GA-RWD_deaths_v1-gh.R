#*[-----------------------------------------------------------------------------------------------]*#
#*[ Objectives : This program applies our changepoint method to COVID-19 daily deaths             ]*#
#*[              using a random walk model with piecewise drifts.                                 ]*#
#*[ Last update: 2021-04-28                                                                       ]*#
#*[ Author     : Jaechoul Lee                                                                     ]*#
#*[-----------------------------------------------------------------------------------------------]*#

### Setup data input and output directories
WD.lib <- c("L:/Home/JaechoulLee/!1Research/Paper/Epidemiology/P02_Covid19/R_library/")
WD.inp <- c("L:/Home/JaechoulLee/!1Research/Paper/Epidemiology/P02_Covid19/Data/")
WD.out <- c("L:/Home/JaechoulLee/!1Research/Paper/Epidemiology/P02_Covid19/Application_R1/Out-deaths_2/")

### Load the R code for the proposed piecewise drift random walk model and GA
source(file=paste(WD.lib,"lib_ga-RWD_v1-2.R",sep=""))

### Required packages
library(tidyverse)
library(doSNOW)
library(foreach)
library(doParallel)

#*[-----------------------------------------------------------------------------------------------]*#
### Data pre-processing
#*[-----------------------------------------------------------------------------------------------]*#

### Read the data trimmed for the study period 2020-03-08 to 2021-02-28
data.STATES <- read.csv(file=paste(WD.inp,"us-states_20210303.csv",sep=""))

names(table(data.STATES$state))

states.m <- c("Alabama","Alaska","Arizona","Arkansas","California",
              "Colorado","Connecticut","Delaware","District of Columbia","Florida",
              "Georgia","Hawaii","Idaho","Illinois","Indiana",
              "Iowa","Kansas","Kentucky","Louisiana","Maine",
              "Maryland","Massachusetts","Michigan","Minnesota","Mississippi",
              "Missouri","Montana","Nebraska","Nevada","New Hampshire",
              "New Jersey","New Mexico","New York","North Carolina","North Dakota",
              "Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island",
              "South Carolina","South Dakota","Tennessee","Texas","Utah",
              "Vermont","Virginia","Washington","West Virginia","Wisconsin",
              "Wyoming")

n.states <- length(states.m)

### Find the first and last recording dates for each state
date.first <- numeric()
date.last  <- numeric()

for (i in 1:n.states) {
  state_i <- states.m[i]                                 # i-th state
  data_i  <- filter(data.STATES,state==state_i)          # i-th state dataset

  n_i <- nrow(data_i)                                    # i-th state sample size
  date.first_i <- data_i$date[1]                         # initial first date
  date.last_i  <- data_i$date[n_i]                       # last date

  date.first[i] <- date.first_i
  date.last[i]  <- date.last_i
}

state.date <- data.frame(state=states.m,date.first,date.last)
state.date

table(state.date$date.first)

### Create dataset with common first and last dates (2020-03-15 to 2021-03-03)
date.first_target <- "2020-03-15"
data.STATES_target <- numeric()

for (i in 1:n.states) {
  state_i <- states.m[i]                                 # i-th state
  data_i  <- filter(data.STATES,state==state_i)          # i-th state dataset

  date.first_i <- data_i$date[1]                         # initial first date

  if (date.first_i <= date.first_target) {
    data_i.target <- filter(data_i,date>=date.first_target)
  } else {
    fips_i  <- data_i$fips[1]
    n.dif_i <- as.numeric(as.Date(date.first_i)-as.Date(date.first_target))

    data_i.early <- data.frame(date=as.character(as.Date(date.first_target)+0:(n.dif_i-1)),
                               state=state_i,
                               fips=fips_i,
                               cases=0,
                               deaths=0)

    data_i.target <- rbind(data_i.early,data_i)
  }

  data.STATES_target <- rbind(data.STATES_target,data_i.target)
}

dim(data.STATES_target)

### Cumulative deaths in the last day
death.last <- numeric()
for (i in 1:n.states) {
  state_i <- states.m[i]                                 # i-th state
  data_i  <- filter(data.STATES_target,state==state_i)   # i-th state dataset

  n_i.all <- nrow(data_i)                                # sample size in daily data
  death.last_i <- filter(data_i,date==data_i$date[n_i.all])

  death.last <- rbind(death.last,death.last_i)
}
death.last

filter(death.last,deaths<=720)                           # only for states <= 720(=2*n) deaths

        date   state fips cases deaths
1 2021-03-03  Alaska    2 58574    291
2 2021-03-03  Hawaii   15 27549    438
3 2021-03-03   Maine   23 45091    705
4 2021-03-03 Vermont   50 15487    207
5 2021-03-03 Wyoming   56 54616    682

#*[-----------------------------------------------------------------------------------------------]*#
### Changepoint estimation via GA method
#*[-----------------------------------------------------------------------------------------------]*#

### Setup parallel backend to use multiple cores
cores <- detectCores()
cl <- makeCluster(cores-1)
registerDoParallel(cl)

### Record run time
print(paste("#==========[  GA began at",Sys.time()," ]==========#"))

### GA for changepoints detection
GA.sol <- foreach (i=1:n.states,.packages='tidyverse') %dopar% {
  state_i <- states.m[i]                                 # i-th state
  data_i  <- filter(data.STATES_target,state==state_i)   # i-th state dataset

  n_i.all <- nrow(data_i)                                # sample size in daily data

  if (data_i$deaths[n_i.all] > 720) {                    # [!CAUTION!] only for states > 720(=2*n) deaths, excluding
                                                         # Alaska,Hawaii,Maine,Vermont,Wyoming
    t.1_i.loc <- 1+3                                     # location of first day for 7-day moving-average
    t.n_i.loc <- n_i.all-3                               # location of last day for 7-day moving-average

    t.1_i <- as.Date(data_i$date[t.1_i.loc])             # first day for 7-day moving-average: 2020-03-18
    t.n_i <- as.Date(data_i$date[t.n_i.loc])             # last day for 7-day moving-average:  2021-02-28

    # Set up the series for 7-day moving-average
    x.t_all <- c(data_i$deaths[1],diff(data_i$deaths))   # daily count of new deaths in data
    x.t_all.mis <- ifelse(x.t_all>=0,x.t_all,NA)         # NA if daily count is < 0
    m.t_all <- mov.avg_q(x.t_all.mis,q=3)                # 7(=2q+1)-day moving-average

    t.t <- as.Date(data_i$date[t.1_i.loc:t.n_i.loc])     # day for study period
    y.t <- data_i$deaths[t.1_i.loc:t.n_i.loc]            # daily count of cumulative deaths
    x.t <- x.t_all.mis[t.1_i.loc:t.n_i.loc]              # daily count of new deaths, with negative as missing
    m.t_0 <- m.t_all[t.1_i.loc:t.n_i.loc]                # 2q+1=7 day moving-average

    m.t <- m.t_0                                         # if m.t_0 < 0.25, m.t is set to be around 0.25
    for (t in 1:length(m.t_0)) {
      if (m.t_0[t] < 0.25) m.t[t] <- rnorm(1,mean=0.25,sd=0.01)
    }
    lm.t <- log(m.t)                                     # use log(m.t) for change in growth rate

    # GA changepoint estimation                          # use log(m.t) for change in growth rate
    ga.out <- ga.cpt_ts(y=lm.t,fitness=pnllik.MDL_RWD,gen.size=200,max.itr=150,p.mut=0.05,
                        seed=10*(i-1)+1000*(i==10)+543,is.graphic=FALSE,is.print=FALSE,is.export=FALSE)

    ga.sol <- ga.out$solution                            # GA estimated changepoints
    ga.mdl <- ga.out$val.sol[length(ga.out$val.sol)]     # optimized value of MDL
    ga.day <- ga.sol[-1]+t.t[1]-1                        # dates of changepoints

    ### ML estimation with GA estimated changepoints
    out.ga <- fit.RWD(lm.t,cp=ga.sol)                    # use log(m.t) for change in growth rate

    delta.ga <- out.ga$delta
    sterr.ga <- out.ga$delta.se
    sigma.ga <- out.ga$sigma
    lnlik.ga <- out.ga$loglik

    ml.out <- cbind(delta.ga,sterr.ga,sigma.ga,lnlik.ga)

    list(state=state_i,ga.sol=ga.sol,ga.mdl=ga.mdl,ga.day=ga.day,ml.est=ml.out)
  }                                                      # end of if-statement
}

### Record run time
print(paste("#==========[  GA ended at",Sys.time()," ]==========#"))

### Print GA changepoint and ML estimationi results
GA.sol

### End parallel backend
stopCluster(cl)

### Save to files
t.n_i <- "2021-02-28"                                    # last day for 7-day moving-average (3-03 minus 3)

capture.output(GA.sol,file=paste(WD.out,"out_ga-RWD_deaths_",t.n_i,".txt",sep=""))
save(GA.sol,file=paste(WD.out,"out_ga-RWD_deaths_",t.n_i,".RData",sep=""))

#*[-----------------------------------------------------------------------------------------------]*#
### GA results summary
#*[-----------------------------------------------------------------------------------------------]*#

### Load GA changepoint analysis results
t.1_i <- "2020-03-18"                                    # first day for 7-day moving-average: 3-(15+3)
t.n_i <- "2021-02-28"                                    # last day for 7-day moving-average: 3-(03-3)

load(file=paste(WD.out,"out_ga-RWD_deaths_",t.n_i,".RData",sep="")) # Loaded is GA.sol
GA.sol

### GA changepoint days
GA.dates <- data.frame()
 
for (i in 1:n.states) {
  ga.out_i <- GA.sol[[i]]                                # GA result for i-th state
  GA.dates_i <- merge(ga.out_i$state,c(t.1_i,as.character(ga.out_i$ga.day))) # state & its changepoint dates
  GA.dates <- rbind(GA.dates,GA.dates_i)
}

colnames(GA.dates) <- c("state","cpt.date")
GA.dates

### RWD estimates with GA changepoints considered
GA.mlest <- data.frame()

for (i in 1:n.states) {
  ga.out_i <- GA.sol[[i]]                                # GA result for i-th state
  GA.mlest_i <- ga.out_i$ml.est                          # RWD estimates
  GA.mlest <- rbind(GA.mlest,GA.mlest_i)
}

GA.mlest

### All GA results
GA.result <- data.frame(GA.dates,GA.mlest)
GA.result

### Export output file in .csv format
write.csv(GA.result,file=paste(WD.out,"out_ga-RWD_deaths_",t.n_i,".csv",sep=""))
