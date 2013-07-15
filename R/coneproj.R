coneA=function(y,amat,w=NULL,lst=TRUE){
	if(!is.numeric(y)||length(y)==0)
		stop("y must be a numeric vector of length >= 1 !")
	if(!is.numeric(amat)||!is.matrix(amat)) 
		stop("amat must be a numeric matrix !")
	if(ncol(amat)!=length(y))
		stop("the column number of amat must equal the length of y !")     
	if(!is.null(w)){
		if(!is.numeric(w))
			stop("w must be a numeric vector !")
		if(any(w<0))
			stop("w must be nonnegative !")
		if(length(w)!=length(y))
			stop("w must have the same length as y !")
		else{
			w=diag(w)
			y=sqrt(w)%*%y
			amat=amat%*%solve(sqrt(w))
 			ans=.Call("coneACpp", y, amat, PACKAGE = "coneproj")
			ans$thetahat=solve(sqrt(w), ans$thetahat)
		}
	}
	else{ 
		ans=.Call("coneACpp", y, amat, PACKAGE = "coneproj")
	}
	rslt=list(df=ans$dim, thetahat=ans$thetahat, steps=ans$nrep, message=NULL, convergence=0)
	upper=length(y)^2
	if(rslt$steps>upper){
		warning ("iterations exceed the maximum number allowed !")
		rslt$convergence = 1
		rslt$message="iterations exceed the maximum number allowed !"}
	if(!lst){
		rslt = list(df=ans$dim, thetahat=ans$thetahat, steps=ans$nrep)    
		rslt}
	else{
		rslt}
}


coneB=function(y,delta,vmat,w=NULL,lst=TRUE){
	if(!is.numeric(y)||length(y)==0)
		stop("y must be a numeric vector of length >= 1 !")
	if(!is.numeric(delta)||!is.matrix(delta)) 
		stop("delta must be a numeric matrix !")      
	if(ncol(delta)!=length(y))
		stop("the column number of delta must equal the length of y !")  
	if(!is.numeric(vmat)||!is.matrix(vmat))
		stop("vmat must be a numeric matrix !")
	if(nrow(vmat)!=length(y))
		stop("the row number of vmat must equal the length of y !")
	if(!is.null(w)){
		if(!is.numeric(w))
			stop("w must be a numeric vector !")
		if(any(w<0))
			stop("w must be a nonnegeative vector !")
		if(length(w)!=length(y))
			stop("w must have the same length as y !")
		else{      
			w=diag(w)
			y=sqrt(w)%*%y
			delta=t(sqrt(w)%*%t(delta))
			vmat=sqrt(w)%*%vmat
			ans=.Call("coneBCpp", y, delta, vmat, PACKAGE = "coneproj")# 
			ans$yhat=solve(sqrt(w), ans$yhat)
		}
	}
	else{
		ans=.Call("coneBCpp", y, delta, vmat, PACKAGE = "coneproj")# 
	}
	rslt=list(df=ans$dim, yhat=ans$yhat, steps=ans$nrep, coefs=ans$coefs, message=NULL, convergence=0)
	upper=length(y)^2
	if(rslt$steps>upper){
		warning ("iterations exceed the maximum number allowed !")
		rslt$convergence=1
		rslt$message="iterations exceed the maximum number allowed !"}
	if(!lst){
		rslt=list(df=ans$dim, yhat=ans$yhat, steps=ans$nrep, coefs=ans$coefs)    
		rslt}
	else{
		rslt
	}
}


qprog=function(q,c,amat,b,lst=TRUE){
	if(!is.vector(c))
		c=as.vector(c)
	if(!is.vector(b))
		b=as.vector(b)
	if(ncol(q)!=length(c))
		stop("the length of c should be equal to the column number of q !")
	if(ncol(q)!=ncol(amat))
		stop("the column number of q should be equal to the column number of amat !")
	if(nrow(amat)!=length(b))
		stop("the row number of amat should be equal to the length of b !")
	else{
		ans=.Call("qprogCpp", q, c, amat, b, PACKAGE = "coneproj")
	}
	rslt=list(df=ans$dim, thetahat=ans$thetahat, steps=ans$nrep, message=NULL, convergence=0)
	upper=length(c)^2
	if(rslt$steps>upper){
		warning ("iterations exceed the maximum number allowed !")
		rslt$convergence=1
		rslt$message="iterations exceed the maximum number allowed !"}
	if(!lst){
		rslt=list(df=ans$dim, thetahat=ans$thetahat, steps=ans$nrep)    
		rslt}
	else{
		rslt
	}
}

qrdecomp=function(xm){
	n=length(xm[,1]);m=length(xm[1,]) 
        sm=1e-8
	nq=1
	lv=0
	qm=xm
	i=0
	while(lv==0&i<=m){
		i=i+1
		lv=sum(xm[,i]^2)
	}
	if(lv==0){stop}
	qm[,nq]=xm[,i]/sqrt(lv)
	for(j in (i+1):m){
		res=xm[,j]
		for(k in 1:nq){
			res=res-sum(qm[,k]*xm[,j])*qm[,k]
		}
		lv=sum(res^2)
		if(lv>sm){
			nq=nq+1
			qm[,nq]=res/sqrt(lv)
		}
	}
	qm=qm[,1:nq]
	ans=new.env()
	ans$qm=qm
	ans$rank=nq
	ans
}


constreg=function(y,xmat,amat,test=FALSE)
{
	n=length(y)
	sm=1e-10
	m=length(xmat)/n 
	pmat=xmat%*%solve(crossprod(xmat), t(xmat))	
	umat=chol(crossprod(xmat))	
	uinv=solve(umat)
	atil=amat%*%uinv
	z=t(uinv)%*%t(xmat)%*%y
	ans=.Call("coneACpp", z, atil, PACKAGE = "coneproj")
	bhat=uinv%*%ans$thetahat	
	fhatc=xmat%*%bhat
	fhatuc=pmat%*%y
	ansq=qrdecomp(t(atil))
	atilq=ansq$qm
	z0=z-atilq%*%t(atilq)%*%z
	yhat0=xmat%*%uinv%*%z0
	sse0=sum((y-yhat0)^2)
	sse1=sum((y-fhatc)^2)
	bval=(sse0-sse1)/sse0
	g=qr(atil)
	dim0=m-g$rank
	if(test){
		nloop=1e+4
		if(bval>sm){
			mdist=0:m*0
             		for(iloop in 1:nloop){
				ys=rnorm(n)
				z=t(uinv)%*%t(xmat)%*%ys
				ans=coneA(z,atil)
				l=ans$df+1
				mdist[l]=mdist[l]+1
			}
			mdist=mdist/nloop
			obs=1:(m+1)
			end=max(obs[mdist!=0])
			pval=sum(mdist[1:(dim0+1)])
			for(j in (dim0+2):(end)){
				alp=(j-dim0-1)/2; bet=(n-j+1)/2 
				addl=pbeta(bval,alp,bet)*mdist[j]
				pval=pval+addl
			}
			pval=1-pval
		}else{pval=1}
		rslt=list(constr.fit=fhatc, unconstr.fit=fhatuc, pval=pval, coefs=bhat)		
		rslt		

	}
	else{
		rslt=list(constr.fit=fhatc, unconstr.fit=fhatuc , coefs=bhat)		
		rslt				

	}
}


shapereg=function(y, t, shape, xmat=NULL, w=NULL, test=FALSE){
        delta=makedelta(t,shape)
	nxmat=xmat
	if(is.null(xmat))
		nxmat=matrix(rep(1,length(y)),nrow=length(y),ncol=1)
	if(!is.null(xmat)){
		if(qr(xmat)$rank != ncol(xmat))
			stop("xmat should be full column rank!")
	}
	if(shape==3|shape==4)
		vmat=cbind(nxmat,t)
	else{vmat=nxmat}
	if(qr(vmat)$rank != ncol(vmat))
		stop("vmat should be full column rank!")
	n=length(y)
	ans=coneB(y,delta,vmat,w=w)
	nd=length(delta)/n
	pv=length(vmat)/n
	coefx=ans$coefs[1:pv]
	bvec=ans$coefs[(pv+1):(pv+nd)]
	vhat=vmat%*%coefx	
	theta=t(delta)%*%bvec
	yhat=theta+vhat 
	sse0=sum((y-vhat)^2) 
	sse1=sum((y-ans$yhat)^2)
	bval=(sse0-sse1)/sse0
	dim0=qr(vmat)$rank
	sm=1e-8
	m=length(delta)/n+length(vmat)/n
	if((n-1.5*ans$df)<=0){sdhat2=sse1}
	else{sdhat2=sse1/(n-1.5*ans$df)}
        se2=solve(crossprod(vmat))*sdhat2
        se.beta=rep(0,pv)
	tstat=rep(0,pv)
	pvals.beta=rep(0,pv)
	for(i in 1:pv){
		se.beta[i]=sqrt(se2[i,i])
		tstat[i]=coefx[i]/se.beta[i]
		pvals.beta[i]=1-pt(tstat[i],1.5*ans$df)
	}
	if(test){
		nloop=1e+4
		mdist=0:m*0
		for(iloop in 1:nloop){
			colvmat=ncol(vmat)
			ys=rep(0,n)
			for(i in 1:colvmat){
				ys=ys+vmat[,i]
			}
			ys=ys+rnorm(n)
	    		ans=coneB(ys,delta,vmat)
	    		l=ans$df+1
	    		mdist[l]=mdist[l]+1
		}
		mdist=mdist/nloop
		obs=1:(m+1)
		end=max(obs[mdist!=0])
		if(bval>sm){
  			 pval=sum(mdist[1:(dim0+1)])
  			 for(j in (dim0+2):(end)){
          			alp=(j-dim0-1)/2; bet=(n-j+1)/2	
				addl=pbeta(bval,alp,bet)*mdist[j]
          			pval=pval+addl
        		 }
    			 pval=1-pval
  		}else{pval=1}
		rslt=list(pval=pval, coefs=coefx, constr.fit=yhat, linear.fit=vhat, se.beta=se.beta, pvals.beta=pvals.beta)	
		rslt	
	}
	else{
		rslt=list(coefs=coefx, constr.fit=yhat, linear.fit=vhat, se.beta=se.beta, pvals.beta=pvals.beta)
	}
}

makedelta=function(x,sh){
	n=length(x)
	xs=sort(x)
	xu=1:n*0
	xu=unique(xs)
	n1=length(xu)
	sm=1e-8
	obs=1:n
	if(n1<n){
		bmat=matrix(0,nrow=n-n1,ncol=n)
		row=0
		for(i in 1:n1){
			cobs=obs[x==xu[i]]
			nr=length(cobs)
			if(nr>1){
				for(j in 2:nr){
					row=row+1
					bmat[row,cobs[1]]=-1;bmat[row,cobs[j]]=1
				}
			}
		}	
	}
	if(sh<3){
		amat=matrix(0,nrow=n1-1,ncol=n)
		for(i in 1:(n1-1)){
			c1=min(obs[abs(x-xu[i])<sm]);c2=min(obs[abs(x-xu[i+1])<sm])
			amat[i,c1]=-1;amat[i,c2]=1
		}
		if(sh==2){amat=-amat}
	}else if(sh==3|sh==4){
		amat=matrix(0,nrow=n1-2,ncol=n)
		for(i in 1:(n1-2)){
			c1=min(obs[x==xu[i]]);c2=min(obs[x==xu[i+1]]);c3=min(obs[x==xu[i+2]])
			amat[i,c1]=xu[i+2]-xu[i+1];amat[i,c2]=xu[i]-xu[i+2];amat[i,c3]=xu[i+1]-xu[i]
		}
		if(sh==4){amat=-amat}
	}else if(sh>4){
		amat=matrix(0,nrow=n1-1,ncol=n)
		for(i in 1:(n1-2)){
			c1=min(obs[x==xu[i]]);c2=min(obs[x==xu[i+1]]);c3=min(obs[x==xu[i+2]])
			amat[i,c1]=xu[i+2]-xu[i+1];amat[i,c2]=xu[i]-xu[i+2];amat[i,c3]=xu[i+1]-xu[i]
		}
		if(sh==5){
			c1=min(obs[x==xu[1]]);c2=min(obs[x==xu[2]])
			amat[n1-1,c1]=-1;amat[n1-1,c2]=1
		}
		if(sh==6){
			c1=min(obs[x==xu[n1]]);c2=min(obs[x==xu[n1-1]])
			amat[n1-1,c1]=-1;amat[n1-1,c2]=1
		}
		if(sh==7){
			amat=-amat
			c1=min(obs[x==xu[n1]]);c2=min(obs[x==xu[n1-1]])
			amat[n1-1,c1]=1;amat[n1-1,c2]=-1
		}
		if(sh==8){
			amat=-amat
			c1=min(obs[x==xu[1]]);c2=min(obs[x==xu[2]])
			amat[n1-1,c1]=1;amat[n1-1,c2]=-1
		}
	}
	if(n1<n){
		wmat=matrix(0,nrow=n,ncol=n1)
		for(i in 1:n1){wmat[abs(x-xu[i])<sm,i]=1}
		atil=amat%*%wmat
		delta=t(wmat%*%t(atil)%*%solve(atil%*%t(atil)))
	}else{delta=t(t(amat)%*%solve(amat%*%t(amat)))} 
        dr=length(delta)/n
	if(sh>2&sh<5){
		pr1=cbind(1:n*0+1,x)
		prmat=pr1%*%solve(crossprod(pr1), t(pr1))
		for(i in 1:dr){delta[i,]=delta[i,]-t(prmat%*%delta[i,])}
	}else{
		for(i in 1:dr){delta[i,]=delta[i,]-mean(delta[i,])}
	}
	for(i in 1:dr){delta[i,]=delta[i,]/sqrt(sum(delta[i,]^2))}
	delta
}




















