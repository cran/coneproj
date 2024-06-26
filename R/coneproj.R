############
#polar cone#
############
coneA <- function(y, amat, w = NULL, face = NULL, msg = TRUE){
  if (!is.numeric(y) | length(y) == 0) {
    stop("y must be a numeric vector of length >= 1 !")
  }
  if (!is.numeric(amat) | !is.matrix(amat)) {
  # if (!is.matrix(amat) & !is(amat, "sparseMatrix")){
     stop("amat must be a numeric matrix !")
  # stop("amat must be a matrix or a sparse matrix from the Matrix package!")
  }
  if (ncol(amat) != length(y)) {
    stop("the column number of amat must equal the length of y !")
  }     
#new:
  if (!is.null(face)) {
     obs <- 1:nrow(amat)
     if (!all(face %in% obs)) {
         warning("Edges of the starting face must be a subset of all edges ! User-defined face is not used !")
         face <- matrix(numeric(0))
     } else {
         if (!is.numeric(face)) {
             stop("face must be numeric !")
         }
         if (is.vector(face)) {
            face <- as.matrix(face, ncol=1)
         }
     }
  }
  if (is.null(face)) {
      face <- matrix(numeric(0))
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
      ans <- .Call("coneACpp", y, amat, face, PACKAGE = "coneproj")
      if (ans$nrep > (length(y)^2 - 1)) {
        if (msg) {
          print (paste("Fail to converge in coneproj!Too many steps! Number of steps:", ans$nrep))
        }
      }
      if (ans$nrep > length(y)^2) {
        print ("Fail to converge in conerpoj!nrep > n^2 !")
      }
      ans$thetahat <- solve(sqrt(w), ans$thetahat)
    }
  } else { 
    ans <- .Call("coneACpp", y, amat, face, PACKAGE = "coneproj")
    if (ans$nrep > (length(y)^2 - 1)) {
      if (msg) {
        print (paste("Fail to converge in coneproj!Too many steps! Number of steps:", ans$nrep))
      }
    }
    if (ans$nrep > length(y)^2) {
      print ("Fail to converge in conerpoj!nrep > n^2 !")
    }
  }
  face <- which(ans$h == 1)
  rslt <- list(df = ans$dim, thetahat = ans$thetahat, steps = ans$nrep, xmat = ans$xmat, face = face, pmat = ans$pmat)
  attr(rslt, "sub") <- "coneA"
  class(rslt) <- "coneproj"
  return (rslt) 
}

#################
#constraint cone#
#################
coneB <- function(y, delta, vmat = NULL, w = NULL, face = NULL, msg = TRUE){
  if (!is.numeric(y) | length(y) == 0) {
    stop("y must be a numeric vector of length >= 1 !")
  }
  if (!is.numeric(delta) | !is.matrix(delta)) { 
    stop("delta must be a numeric matrix !")      
  }
#new:
  if (nrow(delta) != length(y)) {
    stop("the row number of delta must equal the length of y !")
  }
  if (!is.null(vmat)) {
   # if (!is.numeric(vmat) | !is.matrix(vmat)) {
   #   stop("vmat must be a numeric matrix !")
   # }
    if (!is.numeric(vmat)) {
       stop("vmat must be numeric !")
    }
    if (is.vector(vmat)) {
       vmat <- as.matrix(vmat, ncol = 1)
    }
    if (nrow(vmat) != length(y)) {
      stop("the row number of vmat must equal the length of y !")
    }
# if vmat != NULL, we make a copy of vmat 
    nvmat <- vmat
    np <- ncol(vmat)
  }
  if (is.null(vmat)) {
# if vmat == NULL, we make a matrix of a numeric(0) element to replace vmat
    nvmat <- matrix(numeric(0))
    np <- 0
  }
#new:
  if (!is.null(face)) {
     obs <- 1:(ncol(delta)+np)
     if (!all(face %in% obs)) {
         warning("Edges of the starting face must be a subset of all edges ! User-defined face is not used !")
         face <- matrix(numeric(0))
     } else {
         if (!is.numeric(face)) {
             stop("face must be numeric !")
         }
         if (is.vector(face)) {
            face <- as.matrix(face, ncol=1)
         }
     }
  }
  if (is.null(face)) {
      face <- matrix(numeric(0))
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
        delta <- sqrt(w) %*% delta
        if (!is.null(vmat)) {
            nvmat <- sqrt(w) %*% nvmat
        }
        ans <- .Call("coneBCpp", y, t(delta), nvmat, face, PACKAGE = "coneproj")
        if (ans$nrep > (length(y)^2 - 1)) {
          if (msg) {
            print (paste("Fail to converge in coneproj!Too many steps! Number of steps:", ans$nrep))
          }
        }
        ans$yhat <- solve(sqrt(w), ans$yhat)
      }
    } else {
      ans <- .Call("coneBCpp", y, t(delta), nvmat, face, PACKAGE = "coneproj")
      if (ans$nrep > (length(y)^2 - 1)) {
        if (msg) {
          print (paste("Fail to converge in coneproj!Too many steps! Number of steps:", ans$nrep))
        }
      }
    }
  face <- which(ans$coefs > 0L)
  rslt <- list(df = ans$dim, yhat = ans$yhat, steps = ans$nrep, coefs = ans$coefs, face = face)
  attr(rslt, "sub") <- "coneB"
  class(rslt) <- "coneproj"
  return (rslt) 
}

#######################
#quadratic programming#
#######################
qprog <- function(q, c, amat, b, face = NULL, msg = TRUE) {
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
#new:
  if (!is.null(face)) {
     obs <- 1:nrow(amat)
     if (!all(face %in% obs)) {
         warning("Edges of the starting face must be a subset of all edges ! User-defined face is not used !")
         face <- matrix(numeric(0))
     } else {
         if (!is.numeric(face)) {
             stop("face must be numeric !")
         }
         if (is.vector(face)) {
            face <- as.matrix(face, ncol=1)
         }
     }
  }
  if (is.null(face)) {
      face <- matrix(numeric(0))
  }
  if (nrow(amat) != length(b)) {
    stop("the row number of amat should be equal to the length of b !")
  } else {
    ans <- .Call("qprogCpp", q, c, amat, b, face, PACKAGE = "coneproj")
    if (ans$nrep > (length(c)^2 - 1)) {
       if (msg) {
          print (paste("Fail to converge in coneproj!Too many steps! Number of steps:", ans$nrep))
        }
    }
  }
  face <- which(ans$h == 1)
  rslt <- list(df = ans$dim, thetahat = ans$thetahat, steps = ans$nrep, xmat = ans$xmat, face = face)
  attr(rslt, "sub") <- "qprog"
  class(rslt) <- "coneproj"
  return (rslt) 
}

##################################################################
# computes Q of QR decomposition; returns col rank of a matrix   #
##################################################################
qrdecomp <- function(xm) {
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
  return (ans)
}

######################################
#parametrically restricted regression#
######################################
constreg <- function(y, xmat, amat, w = NULL, test = FALSE, nloop = 1e+4) {
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
      }
  }
  #make the projection matrix for the unconstrained alternative
  if (is.null(w)) {
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
  #ansq <- qrdecomp(t(atil))
  #atilq <- ansq$qm
  #z0 <- z - atilq %*% t(atilq) %*% z
#new:
  vmat = qr.Q(qr(t(atil)), complete = TRUE)[, -(1:(qr(t(atil))$rank)), drop = FALSE]
  pvmat = vmat %*% solve(crossprod(vmat), t(vmat))
  z0 = pvmat %*% z 
  yhat0 <- xmat %*% uinv %*% z0
#new: include weighted case: already checked w
  if (!is.null(w)) { 
    sse0 <- sum(w * (y - yhat0)^2)
    sse1 <- sum(w * (y - fhatc)^2)
  } else {
    sse0 <- sum((y - yhat0)^2)
    sse1 <- sum((y - fhatc)^2)
  }
  bval <- (sse0 - sse1) / sse0
  g <- qr(atil)
  dim0 <- m - g$rank
  #-----------------------------------------------------------------
  #new: mixture variance-covariance matrix added in 2023; start
  #-----------------------------------------------------------------
  sighat <- sse1 / (n - ans$df)
  sdhat <- sighat^.5
  k <- length(bhat)
  mcovp <- matrix(0, nrow = k, ncol = k)
  ik <- diag(k)
  nsim <- 1e+3
  for(isim in 1:nsim){
    ysim <- xmat %*% bhat + rnorm(n, 0, sd = sdhat)
    #zsim <- slinv%*%t(umat)%*%t(xmat)%*%ysim/sqrt(n)
    if (!is.null(w)) {
      zsim <- t(uinv) %*% t(xmat) %*% w %*% ysim
    } else {
      zsim <- t(uinv) %*% t(xmat) %*% ysim
    }
    asim <- coneA(zsim, atil)
    if(length(asim$face) == 0){
      mcovp <- mcovp + ik
    } else if(length(asim$face) > 0){
      aj = atil[asim$face, ]
      if(length(aj) == k){aj = matrix(aj, nrow = 1)}
      mcovp = mcovp + ik - t(aj) %*% solve(aj %*% t(aj)) %*% aj
    }		
  }
  mcovp <- mcovp / nsim * sighat
  #mcov <- umat%*%slinv%*%mcovp%*%slinv%*%t(umat)/n
  mcov <- uinv %*% mcovp %*% t(uinv)
  level <- .95
  t.mult <- qt((1 - level)/2, df = (n - ans$df), lower.tail=FALSE)
  pmc <- t.mult * sqrt(diag(mcov))
  uppc.bhat <- bhat + pmc
  lowc.bhat <- bhat - pmc
  
  #I add the following: c.i. for fhat
  hl <- t.mult * sqrt(diag(xmat %*% mcov %*% t(xmat)))
  uppc.fhat <- fhatc + hl
  lowc.fhat <- fhatc - hl
  
  #I add the following: t tests for bhat
  se.beta <- sqrt(diag(mcov))
  tstat <- bhat / se.beta
  if ((n - 1.5*ans$df) <= 0) {
    pvals.beta <- 2 * (1 - pt(abs(tstat),  ans$df))
    warning('Effective degrees of freedom is close to the number of observations! 
            Inference about parametric covariates is not reliable!')
  } else {
    pvals.beta <- 2 * (1 - pt(abs(tstat),  n - 1.5*ans$df))
  }
  #-----------------------------------------------------------------
  #new: mixture variance-covariance matrix added in 2023; end
  #-----------------------------------------------------------------
  #find the p-value for E01 test
  if (test) {
    if (bval > sm) {
      mdist <- 0:m*0
      for (iloop in 1:nloop) {
        ys <- rnorm(n)
#new: include weighted case: already checked w
        if (!is.null(w)) {
	        z <- t(uinv) %*% t(xmat) %*% w %*% ys
	      } else {
          z <- t(uinv) %*% t(xmat) %*% ys
        }
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
    #rslt <- list(constr.fit = fhatc, unconstr.fit = fhatuc, pval = pval, coefs = bhat)
    rslt <- list(constr.fit = fhatc, unconstr.fit = fhatuc, pval = pval, coefs = bhat,
                 uppc.bhat = uppc.bhat, lowc.bhat = lowc.bhat, uppc.fhat = uppc.fhat, 
                 lowc.fhat = lowc.fhat, acov = mcov, se.beta = se.beta, pvals.beta = pvals.beta, tstat = tstat)		
  } else {
    #rslt <- list(constr.fit = fhatc, unconstr.fit = fhatuc , coefs = bhat)
    rslt <- list(constr.fit = fhatc, unconstr.fit = fhatuc , coefs = bhat, 
                 uppc.bhat = uppc.bhat, lowc.bhat = lowc.bhat, uppc.fhat = uppc.fhat, 
                 lowc.fhat = lowc.fhat, acov = mcov, se.beta = se.beta, pvals.beta = pvals.beta, tstat = tstat)			
  }
  attr(rslt, "sub") <- "constreg"
  #class(rslt) <- "coneproj"	
  #new: add in 2023 
  class(rslt) <- c("constreg", "coneproj")
  return (rslt)	
}

###################
#new: add in 2023 #
###################
vcov.constreg <- function(object,...){
  v <- object$acov
  dimnames(v) = list(names(coef(object)), names(coef(object)))
  return(v)
}

############################
#shape-restricted routine  #
#to be consistent with cgam#
############################
shapereg <- function(formula, data = NULL, weights = NULL, test = FALSE, nloop = 1e+4) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  shape <- NULL
  t <- NULL
  tnm <- NULL
  zmat <- NULL; zid <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_param <- NULL; is_fac <- NULL; vals <- NULL; st <- 1; ed <- 1
  for (i in 2: ncol(mf)) {
    if (is.numeric(attributes(mf[,i])$shape)) {
      shape <- c(shape, attributes(mf[,i])$shape)
      t <- cbind(t, mf[,i])
      tnm <- c(tnm, attributes(mf[,i])$nm)
    }
    if (is.null(attributes(mf[,i])$shape)) {
      if (!is.null(names(mf)[i])) {
        znms <- c(znms, names(mf)[i])
      }
      if (!is.matrix(mf[,i])) {
          zid <- c(zid, i)
	  is_param <- c(is_param, TRUE)
          if (is.factor(mf[,i])) {
	    is_fac <- c(is_fac, TRUE)
	    ch_char <- suppressWarnings(is.na(as.numeric(levels(mf[, i]))))
            if (any(ch_char)) {
	      vals <- c(vals, unique(levels(mf[, i]))[-1])
            } else {
	      vals <- c(vals, as.numeric(levels(mf[, i]))[-1])
	    }
            nlvs <- length(attributes(mf[,i])$levels)
	    ed <- st + nlvs - 2 
	    zid1 <- c(zid1, st)
	    zid2 <- c(zid2, ed)
	    st <- st + nlvs - 1
	    zmat0 <- model.matrix(~ mf[, i])[, -1, drop = FALSE]
	    zmat <- cbind(zmat, zmat0)
          } else {
	    is_fac <- c(is_fac, FALSE)
            zmat <- cbind(zmat, mf[, i])
	    ed <- st
            zid1 <- c(zid1, st)
	    zid2 <- c(zid2, ed) 
	    st <- st + 1
	    vals <- c(vals, "")
         }
      } else {
	  is_param <- c(is_param, FALSE)
          is_fac <- c(is_fac, FALSE)
	  zmat0 <- mf[, i]
	  mat_cols <- ncol(zmat0)
	  mat_rm <- NULL
	  #rm_num <- 0
	  for (irm in 1:mat_cols) {
       	  	if (all(round(diff(zmat0[, irm]), 8) == 0)) {
                	mat_rm <- c(mat_rm, irm)
          	}
   	  }
	  if (!is.null(mat_rm)) {
	  	zmat0 <- zmat0[, -mat_rm, drop = FALSE]
		#rm_num <- rm_num + length(mat_rm)
	  }
	  zmat <- cbind(zmat, zmat0)
	  vals <- c(vals, 1)
	  zid <- c(zid, i)
	  nlvs <- ncol(zmat0) + 1
	  ed <- st + nlvs - 2
	  zid1 <- c(zid1, st)
	  zid2 <- c(zid2, ed)
	  st <- st + nlvs - 1
      }
    }  
  }
  if (!test & nloop > 0) {
	nloop <- 0
  	#print ("nloop > 0, test should be TRUE!")
  }
  ans <- shapereg.fit(y, t, shape, xmat = zmat, w = weights, test = test, nloop = nloop)
  rslt <- list(coefs = ans$coefs, constr.fit = ans$constr.fit, linear.fit = ans$linear.fit, se.beta = ans$se.beta, pval = ans$pval, pvals.beta = ans$pvals.beta, test = ans$test, SSE0 = ans$SSE0, SSE1 = ans$SSE1, shape = shape, tms = mt, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2, tnm = tnm,  ynm = ynm, znms = znms, is_param = is_param, is_fac = is_fac, xmat = zmat)
  rslt$call <- match.call()		
  attr(rslt, "sub") <- "shapereg"
  class(rslt) <- "coneproj"
  return (rslt) 
}

##########################
#core routine of shapereg#
##########################
shapereg.fit <- function(y, t, shape, xmat = NULL, w = NULL, test = FALSE, nloop = 1e+4) {
#find delta for the constraint cone 
    delta <- makedelta(t, shape)
    if (is.null(xmat)) {
      nxmat <- NULL
      if (shape == 3 | shape == 4) {
        vmat <- cbind(rep(1, length(y)), t)
      } else {vmat <- matrix(rep(1, length(y)), ncol = 1)}
    }
    if (!is.null(xmat)) {
       if (is.vector(xmat)) {
        nxmat <- matrix(xmat, ncol = 1)
        if (all(round(diff(nxmat), 8) == 0)) {
        #if (all.equal(round(diff(nxmat), 8), matrix(rep(0, nrow(nxmat)-1), ncol = 1))) {
           nxmat <- NULL
        } 
        if (shape == 3 | shape == 4) {
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
          if (all(round(diff(nxmat[, i]), 8) == 0)) {
          #if (all.equal(round(diff(nxmat[, i]), 8), rep(0, mat_rows-1))){
            mat_rm <- c(mat_rm, i)
          }
        }
        if (!is.null(mat_rm)) {
          nxmat <- nxmat[,-mat_rm,drop=FALSE]
        }
        if (shape == 3 | shape == 4) {
          vmat <- cbind(rep(1, length(y)), nxmat, t)
        } else {vmat <- cbind(rep(1, length(y)), nxmat)}
    }
    if (qr(vmat)$rank != ncol(vmat)) {
      stop("vmat should be full column rank!")
    }
  }
  n <- length(y)
  ans <- coneB(y, t(delta), vmat, w)
  nd <- length(delta) / n
  pv <- length(vmat) / n
  #find coefs for vmat and delta
  coefx <- ans$coefs[1:pv]
  if (shape == 3 | shape == 4) {
    coefb <- coefx[-pv]
  } else {coefb <- coefx}
  bvec <- ans$coefs[(pv + 1):(pv + nd)]
#new code to get se(alpha):
  duse <- bvec > 1e-8
  if (sum(duse) >= 1) {
    delsm <- delta[duse, , drop = FALSE]
#debug: include vmat in bmat, not only z's
    #bmat <- cbind(rep(1, n), nxmat, t(delsm))
    bmat <- cbind(vmat, t(delsm))
  } else {
    #bmat <- cbind(rep(1, n), nxmat) 
    bmat <- vmat
  }
  #find H0 fit
  vhat <- vmat %*% coefx	
  theta <- t(delta) %*% bvec
  #find H1 fit
  yhat <- theta + vhat 
#debug: include weighted case in sse0 and sse1
  if (is.null(w)) {
    w <- 1:n*0 + 1
  }
  sse0 <- sum(w * (y - vhat)^2) 
  sse1 <- sum(w * (y - ans$yhat)^2)
  bval <- (sse0 - sse1) / sse0
  dim0 <- qr(vmat)$rank
  sm <- 1e-8
  m <- length(delta) / n + length(vmat) / n
  #find the approximate covariance matrix for beta 
#new: use ans$df is (n - 1.5 * ans$df) <= 0
  #if ((n - 1.5 * ans$df) <= 0) {
  if ((n - dim0 - 1.5 * ans$df) <= 0) {
    sdhat2 <- sse1 / ans$df 
  } else {
    sdhat2 <- sse1 / (n - dim0 - 1.5 * ans$df)
  }
#new code: se(alpha)
  #se2 <- solve(crossprod(vmat)) * sdhat2
  if (is.null(nxmat)) {
    pb <- 1
  } else {
    pb <- 1 + ncol(nxmat)
  }
#debug: include weighted case in se2  
  #se2 <- solve(crossprod(bmat)) * sdhat2
  se2 <- solve(t(bmat) %*% diag(w) %*% bmat) * sdhat2
  se.beta <- rep(0, pb)
  tstat <- rep(0, pb)
  #find the approximate p-values for beta
  pvals.beta <- rep(0, pb)
#new: consistent with cgam
  #if ((n - 1.5 * ans$df) <= 0) {
  #    warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
  #    #pvals.beta <- NULL
  #} else {
  #    for (i in 1:pb) {
  #      se.beta[i] <- sqrt(se2[i, i])
  #      tstat[i] <- coefb[i] / se.beta[i]
  #      pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]), n - 1.5 * ans$df))
  #    }
  #}
  for (i in 1:pb) {
    se.beta[i] <- sqrt(se2[i, i])
    tstat[i] <- coefb[i] / se.beta[i]
    if ((n - dim0 - 1.5 * ans$df) <= 0) {
      pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]), ans$df))
      warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
    } else {
        pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]), n - dim0 - 1.5 * ans$df))
    }
  }
  #find the p-value for E01 test 
  if (test) {
  #find the mixing parameters for the mixture-of-betas distribution for the E01 test statistic
    #nloop <- nloop
    #set.seed(123)
    mdist <- 0:m*0
    for (iloop in 1:nloop) {
      colvmat <- ncol(vmat)
      ys <- rep(0, n)
      for (i in 1:colvmat) {
        ys <- ys + vmat[,i]
      }
      ys <- ys + rnorm(n)
      ans <- coneB(ys, t(delta), vmat, w)
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
      rslt <- list(pval = pval, coefs = coefb, constr.fit = yhat, linear.fit = vhat, se.beta = se.beta, pvals.beta = pvals.beta, shape = shape, test = test, SSE0 = sse0, SSE1 = sse1)
  } else {
      rslt <- list(pval = NULL, coefs = coefb, constr.fit = yhat, linear.fit = vhat, se.beta = se.beta, pvals.beta = pvals.beta, shape = shape, test = test, SSE0 = sse0, SSE1 = sse1)
    }
  #rslt$call <- match.call()		
  #attr(rslt, "sub") <- "shapereg"
  #class(rslt) <- "coneproj"
  return (rslt)
}

###################
#new fitted method#
###################
fitted.coneproj <- function(object,...) {
  sub <- attributes(object)$sub
  if (sub == "coneA" | sub == "qprog") {
    ans <- object$thetahat
  } else if (sub == "coneB") {
    ans <- object$yhat
  } else if (sub == "constreg") {
    #ans <- list(constr.fit = as.vector(object$constr.fit), unconstr.fit = as.vector(object$unconstr.fit))
    ans <- as.vector(object$constr.fit)
  } else if (sub == "shapereg") {
    #ans <- list(constr.fit = object$constr.fit, linear.fit = object$linear.fit)
    ans <- object$constr.fit
  } else {
    ans <- NULL
  }
  ans
}

#################
#new coef method#
#################
coef.coneproj <- function(object,...) {
  ans <- object$coefs
  ans	
}

###################################################
#find delta for a specific predictor x and a shape#
################################################### 
makedelta <- function(x, sh) {
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
  return (delta)
}

############################
#summary                   #
#to be consistent with cgam#
############################
summary.coneproj <- function(object,...) {
  sub <- attributes(object)$sub
  if (sub == "shapereg") {
  	if (!is.null(object$coefs)) {
  		coefs <- object$coefs
  		se <- object$se.beta
  		tval <- coefs / se
  		pvalbeta <- object$pvals.beta
  		pval <- object$pval
  		n <- length(coefs)
  		sse0 <- object$SSE0
  		sse1 <- object$SSE1
  		zid <- object$zid
  		shape <- object$shape
  		#zid1 <- object$zid1 - 1 - length(shape)
  		#zid2 <- object$zid2 - 1 - length(shape)
  		zid1 <- object$zid1
  		zid2 <- object$zid2
  		tms <- object$tms
  		zmat <- object$xmat
  		is_param <- object$is_param
  		is_fac <- object$is_fac
  		vals <- object$vals
  		test <- object$test
  #new:
  		rslt1 <- rslt2 <- NULL
  		if (!is.null(pvalbeta)) {
  			rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
  			rownames(rslt1)[1] <- "(Intercept)"
  			if (n > 1) {
  				lzid <- length(zid1)
  				for (i in 1:lzid) {
  					pos1 <- zid1[i]; pos2 <- zid2[i]
  					for (j in pos1:pos2) {
  						if (!is_param[i]) {
  							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")							
  						} else {
  							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")					
  						}	
  					}
  				}
  			}
  			rslt1 <- as.matrix(rslt1)
  		} 
  		if (!is.null(sse0) & !is.null(sse1)) {
  			#rslt2 <- cbind(SSE.Linear = sse0, SSE.Full = sse1)
  #new:
  			rslt2 <- data.frame("SSE.Linear" = sse0, "SSE.Full" = sse1)
  			if (test) {
  			      PVAL <- object$pval
  			      rslt2 <- cbind(rslt2, "P.value.f(t)" = PVAL)
  		        } 
  			rownames(rslt2)[1] <- ""
  #new:
  			if (!is.null(rslt1)) {
  				ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs) 
  			} else {
  				ans <- list(call = object$call, residuals = rslt2, zcoefs = coefs) 
  			}
  			attr(ans, "sub") <- sub
  			if (!is.null(rslt1)) {
  				class(ans) <- "summary.coneproj"
  			}
  			ans
  		} else {
  			if (!is.null(rslt1)) {
  				ans <- list(call = object$call, coefficients = rslt1, zcoefs = coefs)
  			} else {
  				ans <- list(call = object$call, zcoefs = coefs)
  			}
  			attr(ans, "sub") <- sub
  			if (!is.null(rslt1)) {
  				class(ans) <- "summary.coneproj"
  			}			
  			ans
  		}
  	} else {
  		ans <- list(zcoefs = object$coefs)
  		attr(ans, "sub") <- sub
  		class(ans) <- "summary.coneproj"
  		ans
  	}
  } else if (sub == "constreg") {
    #new: added in 2023
    coefs <- object$coefs
    se <- object$se.beta
    tval <- object$tstat 
    pvalbeta <- object$pvals.beta
    rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
    rslt1 <- as.matrix(rslt1)
    colnames(rslt1) <- c("Estimate", "Std.Err", "t value", "Pr(>t)")
    #temp!
    rownames(rslt1) <- paste0("x", 1:nrow(rslt1))
    ans <- list(call = object$call, coefficients = rslt1,  zcoefs = coefs) 
    attr(ans, "sub") <- sub
    class(ans) <- "summary.coneproj"
    ans
  } else {
	  ans <- list(coefs = object$coefs, fhat = fitted(object))
  	return (ans)
  }
}

###############
#print.summary#
###############
print.summary.coneproj <- function(x,...)
{
  sub <- attributes(x)$sub
  if (sub %in% c("shapereg", "constreg")){
    if (sub == "shapereg") {
      cat("Call:\n")
      print(x$call)
      cat("\n")
    }
    cat("Coefficients:")
    cat("\n")
    op <- options(show.coef.Pvalues = FALSE)
    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
    #cat("==============================================================", "\n")
    #cat("Call:\n")
    #print(x$call)
    #cat("\n")
    if (sub == "shapereg") {
      cat("+--------------------------------------+\n")
      cat("|         SSE.Linear vs SSE.Full       |\n")
      cat("+--------------------------------------+\n")
      printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
    }
  } else {
   #print ("not shapereg!") 
    print (x)
  }
}

###########################################
#suppose the columns of a matrix are edges#
#we check the irreducibility of the matrix#
###########################################
#debugged: must set addrownums = FALSE in order to avoid row name errors, see ?tail in R for detail
check_irred <- function(mat) {
  if (!is.matrix(mat)) {
    stop ("only a matrix can be checked irreducibility!")
  }
#new:
  mat <- t(mat)
  n <- nrow(mat)
  sm <- 1e-8
  nmat <- mat
  base <- mat
  hd <- head(mat, 1)
  #tl <- tail(mat, -1, addrownums = FALSE)
  tl <- tail(mat, -1, keepnums = FALSE)
  id <- 1
  #m <- n
  rm_id <- NULL
#new:
  eq_id <- NULL
  eq <- FALSE 
  eq_num <- 0
  while (nrow(tl) >= 1 & nrow(base) >= 1) {
     ans <- coneB(hd, t(tl))
     hd0 <- hd
     if (all(round(as.vector(ans$yhat), 8) == round(as.vector(hd), 8))) {
       hd <- head(tl, 1)
       #tl <- tail(tl, -1, addrownums = FALSE)
       tl <- tail(tl, -1, keepnums = FALSE)
       rm_id <- c(rm_id, id)
     } else {
         ans <- coneB(-hd, t(tl))
         if (all(round(as.vector(ans$yhat), 8) == round(as.vector(-hd), 8))) {
           eq <- TRUE
           eq_num <- eq_num + 1
#new:
	   	     eq_id <- c(eq_id, id)
           hd <- head(tl, 1)
           #tl <- rbind(tail(tl, -1, addrownums = FALSE), hd0)
           tl <- rbind(tail(tl, -1, keepnums = FALSE), hd0)
         } else {
           hd <- head(tl, 1)
           #tl <- rbind(tail(tl, -1, addrownums = FALSE), hd0)   
           tl <- rbind(tail(tl, -1, keepnums = FALSE), hd0)
         }
     } 
     id <- id + 1
     #m <- m - 1
     base <- base[-1, ,drop = FALSE]
  }
  if (!is.null(rm_id)) {
    nmat <- nmat[-rm_id, ,drop = FALSE]
  }
  if (!eq && is.null(rm_id)) {
     print ("edges are irreducible!")
  }
#debug: should devide eq_num by 2, count twice
  #eq_num <- eq_num / 2
  rslt <- list(edge = nmat, reducible = rm_id, equal = eq_id)
  return (rslt)
}
 
#######################
#eight shape functions#
####################### 
incr <- function(x, numknots = 0, knots = 0, space = "E") 
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 1
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

decr <- function(x, numknots = 0, knots = 0, space = "E") 
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 2
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x 
} 

conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 3
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 4
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

incr.conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 5
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

decr.conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 6
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

incr.conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 7
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

decr.conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 8
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}












