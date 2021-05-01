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
fit.RWD <- function(y,cp) {
  n <- length(y)
  n.rgm <- cp[1]+1
  tau <- c(1,cp[-1],n+1)

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
penalty.MDL <- function(y,cp) {
  n <- length(y)
  m <- cp[1]

  if (m == 0) {
    pnt <- 0
  } else {
    tau.ext <- c(cp[-1],n+1)
    n.r <- numeric(length=m)
    for (i in 1:m) {
      n.r[i] <- sum(!is.na(y[tau.ext[i]:(tau.ext[i+1]-1)]))
    }
    pnt <- log(m+1)+0.5*sum(log(n.r))+sum(log(tau.ext[-1]))
  }

  return(pnt)
}

### Compute the penalized log-likelihood with MDL penalty
pnllik.MDL_RWD <- function(y,cp) {
  pnllik.MDL <- -2*sum(fit.RWD(y=y,cp=cp)$loglik)+penalty.MDL(y=y,cp=cp)
  return(pnllik.MDL)
}

### GA for time series
ga.cpt_ts <- function(y,fitness,gen.size,max.itr,p.mut,seed,is.graphic,is.print,is.export) {
  n <- length(y)

  m.max <- 5
  t.rng <- 10:(n-7)
  Confg <- list()
  Confg.sol <- list()
  Confg.ALL <- list()

  Pnlik <- matrix(0,nrow=max.itr,ncol=gen.size)
  Pnlik.sol <- numeric(length=max.itr)

  if (is.graphic) {
    dev.new(width=12,height=6)
  }

  set.seed(seed)
  for (g in 1:1) {
    if (is.print) print(paste("#----------[  Generation =",g,"has begun at",Sys.time()," ]----------#"))

    Confg[[1]] <- as.integer(0)
    j <- 2
    for (k in 1:(gen.size*1000)) {
      m <- rbinom(1,size=m.max,prob=0.4)
      tau <- sort(sample(t.rng,size=m,replace=FALSE))
      chrom <- c(m,tau)
      Confg[[j]] <- as.integer(chrom)

      is.pass <- FALSE
      if (m == 0) {
        is.pass <- TRUE
      } else {
        if (all(diff(c(1,tau,n+1)) > 6)) {
          is.pass <- TRUE
        }
      }

      if (length(unique(Confg[1:j])) == j & is.pass == TRUE) {
        j <- j+1
      }

      if (j > gen.size) break
    }

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
  }

  for (g in 2:max.itr) {
    if (is.print) print(paste("#----------[  Generation =",g,"has begun at",Sys.time()," ]----------#"))

    gen.rank <- rank(-Pnlik[g-1,])
    gen.rank.sum <- sum(gen.rank)

    Confg.pre <- Confg.ALL[[g-1]]
    Confg[[1]] <- Confg.sol[[g-1]]
    Pnlik[g,1] <- Pnlik.sol[g-1]

    j <- 2
    for (k in 2:(gen.size*1000)) {
      loc.prt <- sample(1:gen.size,size=2,replace=FALSE,prob=gen.rank/gen.rank.sum)
      loc.dad <- loc.prt[1]
      loc.mom <- loc.prt[2]
      chrom.dad <- Confg.pre[[loc.dad]]
      chrom.mom <- Confg.pre[[loc.mom]]

      tau_S1 <- sort(union(chrom.dad[-1],chrom.mom[-1]))
      m_S1 <- length(tau_S1)
      if (m_S1 == 0) {
        m_S4 <- rbinom(1,size=2,prob=p.mut)
        tau_S4 <- sort(sample(t.rng,size=m_S4,replace=FALSE))
      } else {
        ran.val_S2 <- runif(m_S1,min=0,max=1)
        tau_S2 <- tau_S1[ran.val_S2 <= 0.5]
        m_S2 <- length(tau_S2)

        ran.val_S3 <- sample(c(-1,0,1),size=m_S2,replace=TRUE,prob=c(0.3,0.4,0.3))
        tau_S3.tmp <- sort(unique(tau_S2+ran.val_S3))
        tau_S3 <- tau_S3.tmp[tau_S3.tmp %in% t.rng]
        m_S3 <- length(tau_S3)

        m_S4.mut <- rbinom(1,size=2,prob=p.mut)
        tau_S4.mut <- sort(sample(t.rng,size=m_S4.mut,replace=FALSE))
        tau_S4 <- sort(unique(c(tau_S3,tau_S4.mut)))
        m_S4 <- length(tau_S4)
      }

      m <- m_S4
      tau <- tau_S4
      chrom <- c(m,tau)

      Confg[[j]] <- as.integer(chrom)

      is.pass <- FALSE
      if (m == 0) {
        is.pass <- TRUE
      } else {
        if (all(diff(c(1,tau,n+1)) > 6)) {
          is.pass <- TRUE
        }
      }

      if (length(unique(Confg[1:j])) == j & is.pass == TRUE) {
        j <- j+1
      }

      if (j > gen.size) break
    }

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
  }

  list(gen.all=Confg.ALL,gen.sol=Confg.sol,val.all=Pnlik,val.sol=Pnlik.sol,solution=chrom.sol)
}
