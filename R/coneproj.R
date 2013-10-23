coneA <- function(y, amat, w = NULL){
  if (!is.numeric(y) || length(y) == 0) {
    stop("y must be a numeric vector of length >= 1 !")
  }
  if (!is.numeric(amat) || !is.matrix(amat)) {
    stop("amat must be a numeric matrix !")
  }
  if (ncol(amat) != length(y)) {
    stop("the column number of amat must equal the length of y !")
  }     
  if (!is.null(w)) {
    if (!is.numeric(w)) {
      stop("w must be a numeric vector !")
    }
    if (any(w < 0)) {
      stop("w must be nonnegative !")
    }
    if (length(w) != length(y)) {
      stop("w must have the same length as y !")
    } else {
      w <- diag(w)
      y <- sqrt(w) %*% y
      amat <- amat %*% solve(sqrt(w))
      ans <- .Call("coneACpp", y, amat, PACKAGE = "coneproj")
      if (ans$nrep > length(y)^2) {
        print ("Fail to converge in conerpoj!nrep > n^2 !")
      }
      ans$thetahat <- solve(sqrt(w), ans$thetahat)
    }
  } else { 
    ans <- .Call("coneACpp", y, amat, PACKAGE = "coneproj")
    if (ans$nrep > length(y)^2) {
      print ("Fail to converge in conerpoj!nrep > n^2 !")
    }
    }
    rslt <- list(df = ans$dim, thetahat = ans$thetahat, steps = ans$nrep)
	#rslt <- list(df = ans$dim, thetahat = ans$thetahat, steps = ans$nrep, message = NULL, convergence = 0)
	#upper <- length(y)^2
	#if(rslt$steps > upper){
	#	warning ("iterations exceed the maximum number allowed !")
	#	rslt$convergence <- 1
	#	rslt$message <- "iterations exceed the maximum number allowed !"}
	#if(!lst){
	#	rslt <- list(df = ans$dim, thetahat = ans$thetahat, steps = ans$nrep)    
	#	rslt}
	#else{
	#	rslt}
}

coneB <- function(y, delta, vmat = NULL, w = NULL){
  if (!is.numeric(y) || length(y) == 0) {
    stop("y must be a numeric vector of length >= 1 !")
  }
  if (!is.numeric(delta) || !is.matrix(delta)) { 
    stop("delta must be a numeric matrix !")      
  }
  if (ncol(delta) != length(y)) {
    stop("the column number of delta must equal the length of y !")  
  }
  if (!is.null(vmat)) {
    if (!is.numeric(vmat) || !is.matrix(vmat)) {
      stop("vmat must be a numeric matrix !")
    }
    if (nrow(vmat) != length(y)) {
      stop("the row number of vmat must equal the length of y !")
  }
# if vmat != NULL, we make a copy of vmat 
    nvmat <- vmat
 }
 if (is.null(vmat)) {
# if vmat == NULL, we make a matrix of a numeric(0) element to replace vmat
   nvmat <- matrix(numeric(0))
 }
  if (!is.null(w)) {
    if (!is.numeric(w)) {
      stop("w must be a numeric vector !")
    }
    if (any(w < 0)) {
      stop("w must be a nonnegeative vector !")
    }
    if (length(w) != length(y)) {
      stop("w must have the same length as y !")
    } else {      
        w <- diag(w)
        y <- sqrt(w) %*% y
	delta <- t(sqrt(w) %*% t(delta))
        if (!is.null(vmat)) {
	  nvmat <- sqrt(w) %*% nvmat
        }
	ans <- .Call("coneBCpp", y, delta, nvmat, PACKAGE = "coneproj")
	ans$yhat <- solve(sqrt(w), ans$yhat)
      }
    } else {
      ans <- .Call("coneBCpp", y, delta, nvmat, PACKAGE = "coneproj")
      }
  rslt <- list(df = ans$dim, yhat = ans$yhat, steps = ans$nrep, coefs = ans$coefs)
	#rslt <- list(df = ans$dim, yhat = ans$yhat, steps = ans$nrep, coefs = ans$coefs, message = NULL, convergence = 0)
	#upper <- length(y)^2
	#if(rslt$steps > upper){
	#	warning ("iterations exceed the maximum number allowed !")
       	#	rslt$convergence <- 1
	#	rslt$message <- "iterations exceed the maximum number allowed !"}
	#if(!lst){
	#	rslt <- list(df = ans$dim, yhat = ans$yhat, steps = ans$nrep, coefs = ans$coefs)    
       	#	rslt}
	#else{
	#	rslt
	#}
}

qprog <- function(q, c, amat, b){
  if (!is.vector(c)) {
    c <- as.vector(c)
  }
  if (!is.vector(b)) {
    b <- as.vector(b)
  }
  if (ncol(q) != length(c)) {
    stop("the length of c should be equal to the column number of q !")
  }
  if (ncol(q) != ncol(amat)) {
    stop("the column number of q should be equal to the column number of amat !")
  }
  if (nrow(amat) != length(b)) {
    stop("the row number of amat should be equal to the length of b !")
  } else {
  ans <- .Call("qprogCpp", q, c, amat, b, PACKAGE = "coneproj")
  }
  rslt <- list(df = ans$dim, thetahat = ans$thetahat, steps = ans$nrep)
	#rslt <- list(df = ans$dim, thetahat = ans$thetahat, steps = ans$nrep, message = NULL, convergence = 0)
	#upper <- length(c)^2
	#if(rslt$steps > upper){
	#	warning ("iterations exceed the maximum number allowed !")
	#	rslt$convergence = 1
	#	rslt$message = "iterations exceed the maximum number allowed !"}
	#if(!lst){
	#	rslt = list(df = ans$dim, thetahat = ans$thetahat, steps = ans$nrep)    
       	#	rslt}
	#else{
	#	rslt
	#}
}

##################################################################
# computes Q of QR decomposition; returns col rank of a matrix   #
##################################################################
qrdecomp <- function(xm){
  n <- length(xm[, 1]); m <- length(xm[1, ]) 
  sm <- 1e-8
  nq <- 1
  lv <- 0
  qm <- xm
  i <- 0
  while (lv == 0 & i <= m) {
    i <- i + 1
    lv <- sum(xm[, i]^2)
  }
  if (lv == 0) {stop}
    qm[, nq] <- xm[, i] / sqrt(lv)
  for (j in (i + 1):m) {
    res <- xm[, j]
    for (k in 1:nq) {
      res <- res - sum(qm[, k] * xm[, j]) * qm[, k]
    }
    lv <- sum(res^2)
    if (lv > sm) {
      nq <- nq + 1
      qm[, nq] <- res / sqrt(lv)
    }
  }
  qm <- qm[, 1:nq]
  ans <- new.env()
  ans$qm <- qm
  ans$rank <- nq
  ans
}

constreg <- function(y, xmat, amat, w = NULL, test = FALSE){
  n <- length(y)
  sm <- 1e-10
  m <- length(xmat) / n 
  if (!is.null(w)) {
    if (!is.numeric(w)) {
      stop("w must be a numeric vector !")
    }
    if (any(w < 0)) {
      stop("w must be nonnegative !")
    }
    if (length(w) != length(y)) {
      stop("w must have the same length as y !")
    } else {
      w <- diag(w)
      #make the projection matrix for the unconstrained alternative
      pmat <- xmat %*% solve(t(xmat) %*% w %*% xmat, t(xmat) %*% w)
      Q <- t(xmat) %*% w %*% xmat
      umat <- chol(Q)
      uinv <- solve(umat)
      atil <- amat %*% uinv
      z <- t(uinv) %*% t(xmat) %*% w %*% y
      ans <- coneA(z, atil)
      #bhat <- uinv %*% ans$thetahat	
      #get constrained fit
      #fhatc <- xmat %*% bhat
      #get unconstrained fit
      #fhatuc <- pmat %*% y
      }
  }
  #make the projection matrix for the unconstrained alternative
  if(is.null(w)){
    pmat <- xmat %*% solve(crossprod(xmat), t(xmat))	
    #use quadratic programming to get bhat
    umat <- chol(crossprod(xmat))	
    uinv <- solve(umat)
    atil <- amat %*% uinv
    z <- t(uinv) %*% t(xmat) %*% y
    ans <- coneA(z, atil)
  }
  bhat <- uinv %*% ans$thetahat	
  #get constrained fit
  fhatc <- xmat %*% bhat
  #get unconstrained fit
  fhatuc <- pmat %*% y
  #use QR decomposition to get fit under H0
  ansq <- qrdecomp(t(atil))
  atilq <- ansq$qm
  z0 <- z - atilq %*% t(atilq) %*% z
  yhat0 <- xmat %*% uinv %*% z0
  sse0 <- sum((y - yhat0)^2)
  sse1 <- sum((y - fhatc)^2)
  bval <- (sse0 - sse1) / sse0
  g <- qr(atil)
  dim0 <- m - g$rank
  #find the p-value for E01 test
  if (test) {
    nloop <- 1e+4
    set.seed(123)
    if (bval > sm) {
      mdist <- 0:m*0
      for (iloop in 1:nloop) {
        ys <- rnorm(n)
        z <- t(uinv) %*% t(xmat) %*% ys
        ans <- coneA(z, atil)
        l <- ans$df + 1
        mdist[l] <- mdist[l] + 1
      }
      mdist <- mdist / nloop
      obs <- 1:(m + 1)
      end <- max(obs[mdist != 0])
      pval <- sum(mdist[1:(dim0 + 1)])
      for (j in (dim0 + 2):(end)) {
        alp <- (j - dim0 - 1) / 2; bet <- (n - j + 1) / 2 
        addl <- pbeta(bval, alp, bet) * mdist[j]
        pval <- pval + addl
      }
      pval <- 1 - pval
    } else {pval <- 1}
    rslt <- list(constr.fit = fhatc, unconstr.fit = fhatuc, pval = pval, coefs = bhat)		
    rslt		
  } else {
    rslt <- list(constr.fit = fhatc, unconstr.fit = fhatuc , coefs = bhat)		
    rslt				
    }
}

shapereg <- function(y, t, shape, xmat = NULL, w = NULL, test = FALSE,...)UseMethod("shapereg")
shapereg <- function(y, t, shape, xmat = NULL, w = NULL, test = FALSE){
  #find delta for the constraint cone 
    delta <- makedelta(t, shape)
    if (is.null(xmat)) {
      if (shape == 3|shape == 4) {
        vmat <- cbind(rep(1, length(y)), t)
      } else {vmat <- matrix(rep(1, length(y)), ncol = 1)}
    }
    if (!is.null(xmat)) {
      if (is.vector(xmat)) {
        nxmat <- matrix(xmat, ncol = 1)
        if (all.equal(diff(nxmat), matrix(rep(0, nrow(nxmat)-1), ncol = 1)) == TRUE) {
           nxmat <- NULL
        } 
        if (shape == 3|shape == 4) {
          vmat <- cbind(rep(1, length(y)), nxmat, t)
        } else {vmat <- cbind(rep(1, length(y)), nxmat)}
      }
      if (is.matrix(xmat)) {
        nxmat <- xmat
        if (qr(nxmat)$rank != ncol(nxmat)) {
          stop("xmat should be full column rank!")
        }
        mat_cols <- ncol(nxmat)
        mat_rows <- nrow(nxmat)
        mat_rm <- NULL
        for (i in 1:mat_cols) {
          if (all.equal(diff(nxmat[, i]), rep(0, mat_rows-1)) == TRUE){
            mat_rm <- c(mat_rm, i)
          }
        }
        if (!is.null(mat_rm)) {
          nxmat <- nxmat[,-mat_rm,drop=FALSE]
        }
        if (shape == 3|shape == 4) {
          vmat <- cbind(rep(1, length(y)) ,nxmat, t)
        } else {vmat <- cbind(rep(1, length(y)), nxmat)}
    }
    if (qr(vmat)$rank != ncol(vmat)) {
      stop("vmat should be full column rank!")
    }
  }
  n <- length(y)
  ans <- coneB(y, delta, vmat, w)
  nd <- length(delta) / n
  pv <- length(vmat) / n
  #find coefs for vmat and delta
  coefx <- ans$coefs[1:pv]
  bvec <- ans$coefs[(pv + 1):(pv + nd)]
  #find H0 fit
  vhat <- vmat %*% coefx	
  theta <- t(delta) %*% bvec
  #find H1 fit
  yhat <- theta + vhat 
  sse0 <- sum((y - vhat)^2) 
  sse1 <- sum((y - ans$yhat)^2)
  bval <- (sse0 - sse1) / sse0
  dim0 <- qr(vmat)$rank
  sm <- 1e-8
  m <- length(delta) / n + length(vmat) / n
  #find the approximate covariance matrix for beta 
  if ((n - 1.5 * ans$df) <= 0) {
    sdhat2 <- sse1
  } else {
    sdhat2 <- sse1 / (n - 1.5 * ans$df)
    }
  se2 <- solve(crossprod(vmat)) * sdhat2
  se.beta <- rep(0, pv)
  tstat <- rep(0, pv)
  #find the approximate p-values for beta
  pvals.beta <- rep(0, pv)
  for (i in 1:pv) {
    se.beta[i] <- sqrt(se2[i, i])
    tstat[i] <- coefx[i] / se.beta[i]
    pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]), n - 1.5 * ans$df))
  }
  #find the p-value for E01 test 
  if (test) {
  #find the mixing parameters for the mixture-of-betas distribution for the E01 test statistic
    nloop <- 1e+4
    set.seed(123)
    mdist <- 0:m*0
    for (iloop in 1:nloop) {
      colvmat <- ncol(vmat)
      ys <- rep(0, n)
      for (i in 1:colvmat) {
        ys <- ys + vmat[,i]
      }
      ys <- ys + rnorm(n)
      ans <- coneB(ys, delta, vmat)
      l <- ans$df + 1
      mdist[l] <- mdist[l] + 1
    }
    mdist <- mdist / nloop
    obs <- 1:(m + 1)
    end <- max(obs[mdist != 0])
    if (bval > sm) {
      pval <- sum(mdist[1:(dim0 + 1)])
      for (j in (dim0 + 2):(end)) {
        alp <- (j - dim0 - 1) / 2; bet <- (n - j + 1) / 2	
        addl <- pbeta(bval, alp, bet) * mdist[j]
        pval <- pval + addl
      }
      pval <- 1 - pval
    } else {pval <- 1}
      rslt <- list(pval = pval, coefs = coefx, constr.fit = yhat, linear.fit = vhat, se.beta = se.beta, pvals.beta = pvals.beta, shape = shape, test = test, SSE0 = sse0, SSE1 = sse1)
      rslt$call <- match.call()	
      class(rslt) <- "shapereg"	
      rslt
  } else {
      rslt <- list(coefs = coefx, constr.fit = yhat, linear.fit = vhat, se.beta = se.beta, pvals.beta = pvals.beta, shape = shape, test = test, SSE0 = sse0, SSE1 = sse1)
      rslt$call <- match.call()	
      class(rslt) <- "shapereg"
      rslt   
    }
}

#find delta for a specific predictor x and a shape 
makedelta <- function(x, sh){
  n <- length(x)
  xs <- sort(x)
  xu <- 1:n*0
  xu <- unique(xs)
  n1 <- length(xu)
  sm <- 1e-8
  obs <- 1:n
  if (n1 < n) {
    bmat <- matrix(0, nrow = n-n1, ncol = n)
    row <- 0
    for (i in 1:n1) {
      cobs <- obs[x == xu[i]]
      nr <- length(cobs)
      if (nr > 1) {
        for (j in 2:nr) {
          row <- row + 1
          bmat[row, cobs[1]] <- -1; bmat[row, cobs[j]] <- 1
        }
      }
    }	
  }
  if (sh < 3) {
    amat <- matrix(0, nrow = n1-1, ncol = n)
    for (i in 1:(n1-1)) {
      c1 <- min(obs[abs(x - xu[i]) < sm]); c2 <- min(obs[abs(x - xu[i + 1]) < sm])
      amat[i, c1] <- -1; amat[i, c2] <- 1
    }
    if (sh == 2) {
      amat <- -amat
    }
  } else if (sh == 3 | sh == 4) {
    amat <- matrix(0, nrow = n1 - 2, ncol = n)
    for (i in 1:(n1-2)) {
      c1 <- min(obs[x == xu[i]]); c2 <- min(obs[x == xu[i + 1]]); c3 <- min(obs[x == xu[i + 2]])
      amat[i, c1] <- xu[i + 2] - xu[i + 1]; amat[i, c2] <- xu[i] - xu[i + 2]; amat[i, c3] <- xu[i + 1] - xu[i]
    }
    if (sh == 4) {
      amat <- -amat
    }
  } else if (sh > 4) {
    amat <- matrix(0, nrow = n1 - 1, ncol = n)
    for (i in 1:(n1 - 2)) {
      c1 <- min(obs[x == xu[i]]); c2 <- min(obs[x == xu[i + 1]]); c3 <- min(obs[x == xu[i + 2]])
      amat[i, c1] <- xu[i + 2] - xu[i + 1]; amat[i, c2] <- xu[i] - xu[i + 2]; amat[i, c3] <- xu[i + 1] - xu[i]
    }
    if (sh == 5) {
      c1 <- min(obs[x == xu[1]]); c2 <- min(obs[x == xu[2]])
      amat[n1 - 1, c1] <- -1; amat[n1 - 1, c2] <- 1
    }
    if (sh == 6) {
      c1 <- min(obs[x == xu[n1]]); c2 <- min(obs[x == xu[n1 - 1]])
      amat[n1 - 1, c1] <- -1; amat[n1 - 1, c2] <- 1
    }
    if (sh == 7) {
      amat <- -amat
      c1 <- min(obs[x == xu[n1]]); c2 <- min(obs[x == xu[n1 - 1]])
      amat[n1 - 1, c1] <- 1; amat[n1 - 1, c2] <- -1
     }
    if (sh == 8) {
      amat <- -amat
      c1 <- min(obs[x == xu[1]]); c2 <- min(obs[x == xu[2]])
      amat[n1 - 1, c1] <- 1; amat[n1 - 1, c2] <- -1
    }
  }
  if (n1 < n) {
    wmat <- matrix(0, nrow = n, ncol = n1)
    for (i in 1:n1) {
      wmat[abs(x - xu[i]) < sm, i] <- 1
    }
    atil <- amat %*% wmat
    delta <- t(wmat %*% t(atil) %*% solve(atil %*% t(atil)))
    } else {
      delta <- t(t(amat) %*% solve(amat %*% t(amat)))
      } 
      dr <- length(delta) / n
  if (sh > 2 & sh < 5) {
    pr1 <- cbind(1:n * 0 + 1, x)
    prmat <- pr1 %*% solve(crossprod(pr1), t(pr1))
    for (i in 1:dr) {
      delta[i, ] <- delta[i, ] - t(prmat %*% delta[i, ])
    }
  } else {
      for (i in 1:dr) {
        delta[i, ] <- delta[i, ] - mean(delta[i, ])
      }
    }
  for (i in 1:dr) {
    delta[i, ] <- delta[i, ] / sqrt(sum(delta[i, ]^2))
  }
  delta
}

summary.shapereg <- function(object,...)
{
  s <- object$shape
  coefs <- object$coefs
  se <- object$se.beta
  tval <- coefs/se
  pvalbeta <- object$pvals.beta
  n <- length(coefs)
  sse0 <- object$SSE0
  sse1 <- object$SSE1
  test <- object$test
  if (s == 3|s == 4){
    rslt1 <- data.frame(Estimate = coefs[1:(n-1)], StdErr = se[1:(n-1)], t.value = tval[1:(n-1)], p.value = pvalbeta[1:(n-1)])
    rownames(rslt1)[1] <- "intercept"
    if (n > 1){
      num <- 2:(n-1)
      for (i in num){rownames(rslt1)[i] <- paste("xmat[,",i-1, "]", sep = "")}
    }
  } else {
       rslt1 <- data.frame(Estimate = coefs, StdErr = se, t.value = tval, p.value = pvalbeta)
       rownames(rslt1)[1] <- "intercept"
       if (n > 1){
         num <- 2:n
         for (i in num){rownames(rslt1)[i] <- paste("xmat[,",i-1, "]", sep = "")}
       }
    }
  rslt1 <- as.matrix(rslt1)
  rslt2 <- cbind(SSE.Linear = sse0, SSE.Full = sse1)
  if (test) {
    PVAL <- object$pval
    rslt2 <- cbind(rslt2, "P.value.f(t)" = PVAL)
  } 
  ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2)
  class(ans) <- "summary.shapereg"
  ans
}

print.summary.shapereg <- function(x,...)
{
   cat("Call:\n")
   print(x$call)
   cat("\n")
   cat("Coefficients:")
   cat("\n")
   printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
   cat("==============================================================")
   cat("\n")
   cat("Call:\n")
   print(x$call)
   cat("\n")
   printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
}


















