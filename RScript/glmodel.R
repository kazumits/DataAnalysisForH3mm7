# General Linear Model (equiv. to lm())
glmodel <- function(Y,X)
{
    X <- as.matrix(X) # design
    Y <- as.matrix(Y) # observation
    p <- ncol(X)
    N <- nrow(Y)
    Ginv <- solve(t(X)%*%X)
    B <- Ginv%*%t(X)%*%Y # MLE of beta
    # null model : intercept only model
    d0 <- colSums(scale(Y,center=TRUE,scale=FALSE)^2)
    d1 <- colSums((Y-X%*%B)^2) # residual deviance
    # unbiased estimation of var(Y)
    ss <- d1/(N-p)
    fstat <- (d0-d1)/(p-1)/ss
    # F-test against intercept only model
    var.p <- pf(fstat,p-1,N-p,lower.tail = FALSE)
    # variance of beta (inverse Fisher information)
    Bvar <- diag(Ginv)%*%t(ss)
    rownames(Bvar) <- colnames(X)
    # Wald test against b = 0
    wald.p <- 2*pt(abs(B/sqrt(Bvar)),df=N-p,lower.tail = FALSE)
    list(
      beta     = drop(B),
      betavar  = drop(Bvar),
      yvar     = drop(ss),
      Ginv     = Ginv,
      deviance = d1,
      fstat    = fstat,
      var.p    = var.p,
      wald.p   = wald.p
    )
}

# Statistical test between contrasts
compareCells <- function(cell,time,chip,B,V)
{
  chiplab <- apply(expand.grid(time,chip),1,function(x) paste(x,collapse='-'))
  sapply(chiplab,function(x) {
    s1 <- paste(cell[1],x,sep='-') 
    s2 <- paste(cell[2],x,sep='-') 
    # Two contrasts from independent models are independent
    2*pnorm(abs(B[,s1]-B[,s2]),0,sqrt(V[,s1]+V[,s2]),lower.tail = FALSE)  
  })
}


