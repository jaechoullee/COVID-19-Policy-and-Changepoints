### Compute (2q+1)-day moving-average series
mov.avg_q <- function(x,q) {
  n <- length(x)
  x.ma <- numeric(length=n)

  for (t in 1:n) {
    x.ma[t] <- mean(x[max(1,t-q):min(t+q,n)],na.rm=TRUE)
  }

  return(x.ma)  
}

### Compute the maximum likelihood estimation for a RW model with stepwise drifts with changepoints given
fit.RWD <- function(y,cp) {                        # y    : data
                                                   # cp   : changepoint chromosome (m; tau_1,...,tau_m)
  n <- length(y)                                   # n    : sample size
  n.rgm <- cp[1]+1                                 # n.rgm: number of regimes
  tau <- c(1,cp[-1],n+1)                           # tau  : changepoint time configulation including 1 & n+1

  y_diff <- c(y[1],diff(y))

  delta.k <- numeric(length=n.rgm)
  sterr.k <- numeric(length=n.rgm)
  sigma.k <- numeric(length=n.rgm)
  lnlik.k <- numeric(length=n.rgm)

  for (k in 1:n.rgm) {
    y_diff.k <- y_diff[tau[k]:(tau[k+1]-1)]
    out.RWD <- arima(y_diff.k,order=c(0,0,0),include.mean=TRUE,method="CSS-ML")

    delta.k[k] <- as.numeric(out.RWD$coef["intercept"])
    sterr.k[k] <- sqrt(out.RWD$var.coef[1])
    sigma.k[k] <- sqrt(out.RWD$sigma2)
    lnlik.k[k] <- out.RWD$loglik
  }

  list(delta=delta.k,delta.se=sterr.k,sigma=sigma.k,loglik=lnlik.k)
}

### Compute the MDL penalty
penalty.MDL <- function(y,cp) {                    # y   : data
                                                   # cp  : changepoint chromosome (m; tau_1,...,tau_m)
  n <- length(y)                                   # n   : sample size
  m <- cp[1]                                       # m   : number of changepoints

  if (m == 0) {
    pnt <- 0
  } else {
    tau.ext <- c(cp[-1],n+1)                       # tau.ext: changepoints in days (tau_1,...,tau_m,N+1)
    n.r <- numeric(length=m)                       # n.r : no. of observations in each regime
    for (i in 1:m) {
      n.r[i] <- sum(!is.na(y[tau.ext[i]:(tau.ext[i+1]-1)]))
    }
    pnt <- log(m+1)+0.5*sum(log(n.r))+sum(log(tau.ext[-1]))
  }

  return(pnt)                                      # smaller is better
}

### Compute the penalized log-likelihood with MDL penalty
pnllik.MDL_RWD <- function(y,cp) {                 # y   : response
                                                   # cp  : changepoint chromosome (m;xi_1,...,xi_m)
  pnllik.MDL <- -2*sum(fit.RWD(y=y,cp=cp)$loglik)+penalty.MDL(y=y,cp=cp)   # fitness: -2*log(L)+penalty
  return(pnllik.MDL)                               # smaller is better
}

### GA for time series
ga.cpt_ts <- function(y,fitness,gen.size,max.itr,p.mut,
                      seed,is.graphic,is.print,is.export) {
  n <- length(y)                                   # sample size

  # Changepoint configuration
  m.max <- 5                                       # max number of possible changepoints
  t.rng <- 10:(n-7)
  Confg <- list()                                  # changepoint configuration for a generation
  Confg.sol <- list()                              # best changepoint configuration for a generation
  Confg.ALL <- list()                              # changepoint configuration for all generations

  Pnlik <- matrix(0,nrow=max.itr,ncol=gen.size)    # penalized likelihood for all chagenpoint configurations
  Pnlik.sol <- numeric(length=max.itr)             # smallest penalized likelihood value for each generation

  if (is.graphic) {
    dev.new(width=12,height=6)
  }

  # Initial generation
  set.seed(seed)
  for (g in 1:1) {
    if (is.print) print(paste("#----------[  Generation =",g,"has begun at",Sys.time()," ]----------#"))

    Confg[[1]] <- as.integer(0)                    # A chromosome of no changepoints is always considered
    j <- 2                                         # loop index for generation
    for (k in 1:(gen.size*1000)) {                 # This loop still works both when n.cpt>=1 and n.cpt=0
      m <- rbinom(1,size=m.max,prob=0.4)           # [!CAUTION!] Adjust prob=0.4 for other settings
      tau <- sort(sample(t.rng,size=m,replace=FALSE))
      chrom <- c(m,tau)                            # changepoint locations (m; tau_1,...,tau_m)
      Confg[[j]] <- as.integer(chrom)

      is.pass <- FALSE
      if (m == 0) {
        is.pass <- TRUE
      } else {
        if (all(diff(c(1,tau,n+1)) > 6)) {         # allow a changepoint grater than 6 time points
          is.pass <- TRUE
        }
      }

      if (length(unique(Confg[1:j])) == j & is.pass == TRUE) {
        j <- j+1                                   # generation increases when (1) a new child chrom is born
      }                                            #                       and (2) the above condition is met  

      if (j > gen.size) break                      # Produce a generation of gen.size
    }                                              # Ending loop in k

    ### Compute penalized log-likelihood for each chromosome
    for (j in 1:gen.size) {
      chrom <- Confg[[j]]

      if (is.print) print(chrom)

      Pnlik[g,j] <- fitness(y=y,cp=chrom)

      if (is.graphic) {
        plot.ts(y,xlab="Day",ylab="Daily number of confirmed cases",col="gray",
                main=paste("Generation",g,"& Child",j,"( PLKHD =",format(Pnlik[g,j],nsmall=3),")"))
        abline(v=chrom[-1],col="blue",lty=2)
      }
    }

    loc.sol <- which(Pnlik[g,] == min(Pnlik[g,]))
    chrom.sol <- Confg[[loc.sol]]
    Confg.sol[[g]] <- chrom.sol
    Confg.ALL[[g]] <- Confg
    Pnlik.sol[g] <- Pnlik[g,loc.sol]

    if (is.export) {
      capture.output(Confg,file=paste(WD.out,sprintf("GA-Gen_%03d.txt",g),sep=""),append=FALSE)
      write.table(t(format(Pnlik[g,],nsmall=12)),file=paste(WD.out,"GA-Pnlik.csv",sep=""),
                  sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE,append=FALSE)
    }
  }                                                # Ending loop in g

  # Next generations from 2 to gen.size
  for (g in 2:max.itr) {
    if (is.print) print(paste("#----------[  Generation =",g,"has begun at",Sys.time()," ]----------#"))

    # Rank chromosomes in the (g-1)th generation
    gen.rank <- rank(-Pnlik[g-1,])
    gen.rank.sum <- sum(gen.rank)

    # Generate g-th generation: the fittest chromosome carries over to next generation
    Confg.pre <- Confg.ALL[[g-1]]
    Confg[[1]] <- Confg.sol[[g-1]]
    Pnlik[g,1] <- Pnlik.sol[g-1]

    j <- 2                                         # index for child in a generation
    for (k in 2:(gen.size*1000)) {
      # Select father and mother chromosomes
      loc.prt <- sample(1:gen.size,size=2,replace=FALSE,prob=gen.rank/gen.rank.sum)
      loc.dad <- loc.prt[1]
      loc.mom <- loc.prt[2]
      chrom.dad <- Confg.pre[[loc.dad]]
      chrom.mom <- Confg.pre[[loc.mom]]

      # Producing child chromosomes
      # Step 1: Combining
      tau_S1 <- sort(union(chrom.dad[-1],chrom.mom[-1]))  # Do not allow identical chagepoint times
      m_S1 <- length(tau_S1)
      if (m_S1 == 0) {
        # Step 2: Thinning (SKIP!!!)
        # Step 3: Shifting (SKIP!!!)
        # Step 4: Mutation
        m_S4 <- rbinom(1,size=2,prob=p.mut)               # [!CAUTION!] Adjust p.mut for other settings
        tau_S4 <- sort(sample(t.rng,size=m_S4,replace=FALSE))
      } else {
        # Step 2: Thinning
        ran.val_S2 <- runif(m_S1,min=0,max=1)
        tau_S2 <- tau_S1[ran.val_S2 <= 0.5]
        m_S2 <- length(tau_S2)

        # Step 3: Shifting
        ran.val_S3 <- sample(c(-1,0,1),size=m_S2,replace=TRUE,prob=c(0.3,0.4,0.3))
        tau_S3.tmp <- sort(unique(tau_S2+ran.val_S3))
        tau_S3 <- tau_S3.tmp[tau_S3.tmp %in% t.rng]       # Changepoints must occur in t.rng
        m_S3 <- length(tau_S3)

        # Step 4: Mutation
        m_S4.mut <- rbinom(1,size=2,prob=p.mut)           # [!CAUTION!] Adjust p.mut for other settings
        tau_S4.mut <- sort(sample(t.rng,size=m_S4.mut,replace=FALSE))
        tau_S4 <- sort(unique(c(tau_S3,tau_S4.mut)))
        m_S4 <- length(tau_S4)
      }

      m <- m_S4                                    # number of changepoints
      tau <- tau_S4
      chrom <- c(m,tau)                            # changepoint locations (m; tau_1,...,tau_m)

      Confg[[j]] <- as.integer(chrom)

      is.pass <- FALSE
      if (m == 0) {
        is.pass <- TRUE
      } else {
        if (all(diff(c(1,tau,n+1)) > 6)) {         # allow a changepoint grater than 6 time points
          is.pass <- TRUE
        }
      }

      if (length(unique(Confg[1:j])) == j & is.pass == TRUE) {
        j <- j+1                                   # generation increases when (1) a new child chrom is born
      }                                            #                       and (2) the above condition is met  

      if (j > gen.size) break                      # Produce a generation of gen.size
    }                                              # Ending loop in k

    # Compute penalized log-likelihood for each chromosome
    for (j in 1:gen.size) {
      chrom <- Confg[[j]]

      if (is.print) print(chrom)

      Pnlik[g,j] <- fitness(y=y,cp=chrom)

      if (is.graphic) {
        plot.ts(y,xlab="Day",ylab="Daily number of confirmed cases",col="gray",
                main=paste("Solution in Generation",g-1,
                           "( PLKHD =",format(Pnlik.sol[g-1],nsmall=3),") vs",
                           "Generation",g,"& Child",j,
                           "( PLKHD =",format(Pnlik[g,j],nsmall=3),")"))
        abline(v=chrom.sol[-1],col="red",lty=1)
        abline(v=chrom[-1],col="blue",lty=2)
      }
    }

    loc.sol <- which(Pnlik[g,] == min(Pnlik[g,]))

    chrom.sol <- Confg[[loc.sol]]
    Confg.sol[[g]] <- chrom.sol
    Confg.ALL[[g]] <- Confg
    Pnlik.sol[g] <- Pnlik[g,loc.sol]

    if (is.print) {
      print(c(k,j))
      print(chrom.sol)
      print(paste("MDL =",format(Pnlik.sol[g],nsmall=3)),quote=FALSE)
    }

    if (is.export) {
      capture.output(Confg,file=paste(WD.out,sprintf("GA-Gen_%03d.txt",g),sep=""),append=FALSE)
      write.table(t(format(Pnlik[g,],nsmall=12)),file=paste(WD.out,"GA-Pnlik.csv",sep=""),
                  sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
    }
  }                                                # Ending loop in g

  list(gen.all=Confg.ALL,gen.sol=Confg.sol,val.all=Pnlik,val.sol=Pnlik.sol,solution=chrom.sol)
}                                                  # Ending function: ga.cpt_ts.Pois

