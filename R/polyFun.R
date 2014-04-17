stepPoly <- function(fullData,deltaStar,tau,lambda,tol){

	obsData <- fullData[which(fullData$delta0*fullData$delta1==1),]
	obsData <- obsData[obsData$x!=0,]
	#obsData <- obsData[which(obsData$x+deltaStar<=obsData$v),]
	
	w <- fullData$w; v <- obsData$v
	x <- obsData$x; n <- nrow(obsData)
	
	function(T1,T2,Q1,xGrid,vGrid,verbose=FALSE,display=FALSE){
		
		if(display){
			dev.set(which=4)
			plot(0,0,type="n",xlim=c(-0.3,0.3),ylim=c(-0.3,0.3),xlab=expression(D[1]/sigma[1]),ylab=expression(D[2]/sigma[2]))
			abline(h=0,col="blue")
			abline(v=0,col="blue")
			#require("ellipse")
			}
		
		lambdaGrid <- lambda(xGrid)
		nx <- length(xGrid); nv <- length(vGrid)
		
		function(Q2){
			#browser()
			delta<- deltaStar
			storage.mode(lambdaGrid)<-"double"
			
			logL <- 0.0
			
			degree <- 1L
			normGrad <- rep(0.0,degree+1)
			epsOpt <- rep(0.0,degree+1)
			
			#z <- .C('tmleStep1',a_delta=delta,a_tau=tau,a_n=n,a_xx=x,a_vv=v,a_lambda=lambdaGrid, a_xGrid=xGrid,a_vGrid=vGrid,a_nx=nx,a_nv=nv,a_T1=T1,a_T2=T2,a_Q2=Q2, a_normGrad=normGrad,a_logL=logL,a_epsOpt=epsOpt)
			
			z <- .C('tmleStep',a_delta=delta,a_tau=tau,a_n=n,a_xx=x,a_vv=v,a_lambda=lambdaGrid, a_xGrid=xGrid,a_vGrid=vGrid,a_nx=nx,a_nv=nv,a_T1=T1,a_T2=T2,a_Q2=Q2, a_normGrad=normGrad,a_logL=logL,a_epsOpt=epsOpt,a_degree=degree)
			
			Q2 <- z$a_Q2
			normGrad <- z$a_normGrad
			logL <- z$a_logL
			epsOpt <- z$a_epsOpt
			
			if(any(Q2<0))stop("Negative Elements in Q2...\a")
			
			if(display){
				dev.set(which=4)
				points(normGrad[1], normGrad[2],col="red")
				#sig <- diag(sqrt(empVar))
				#lines(ellipse(sig, scale = sqrt(empVar)/empMean, centre = empMean, level = 0.95),col="red")#level = 0.001
				#lines(ellipse(sig, centre = empMean, level = 0.95),col="red")
				}

			if(verbose) cat("D: ", normGrad[1]," ", normGrad[2]," |  logLik =",logL," ")
			
			cond <- (max(abs(normGrad))<tol)
			
			return(list("Q2"=Q2,"cond"=cond,"normGrad"= normGrad,"logL"=logL,"epsOpt"=epsOpt))
			
			}
		
		}
	}

polyPsi <- function(fullData, deltaStar,tau){
	
	w <- fullData$w; v <- fullData$v
	
	obsData <- fullData[which(fullData$delta0*fullData$delta1==1),]
	obsData <- obsData[obsData$x!=0,]
	# obsData <- obsData[which(obsData$x+deltaStar<=obsData$v),]
	
	function(T1,T2,Q1,Q2,lambda,xGrid,vGrid){
		
		Sn <- Q1$Sn; Qc <- Q1$Qc; dGn <- Q1$dGn
		
		
		lambdax <- function(x) x*lambda(x);	lambdax2 <- function(x) x^2*lambda(x)
		
		i0lambda<-integrate(lambda,0,tau)$value
		i1lambda<-integrate(lambdax,0,tau)$value
		i2lambda<-integrate(lambdax2,0,tau)$value
		Lambda<-matrix(c(i0lambda,i1lambda,i1lambda,i2lambda),2,2)
		
		lambdaGrid <- lambda(xGrid)
		
		intSn <- function(x){
			t <- c(x,sort(v)[sort(v)>x])
			sum(diff(t)*Sn(t[2:length(t)]))
			}
		intSn <- Vectorize(intSn)
		R1<-intSn(xGrid)/Sn(xGrid)
		R1[is.nan(R1)]<-1

		
		
		R2 <- apply(T1*Q2,2,.intFun,vGrid)/apply(T2*Q2,2,.intFun,vGrid)
		
		R2[is.nan(R2)] <- 1
		R <- log(R1)-log(R2)
		
		R[is.infinite(R)]<-0

		intR <- c(.intFun(R*lambdaGrid,xGrid),.intFun(xGrid*R*lambdaGrid,xGrid))
		
		theta <- solve(Lambda,intR)

		return(theta)
		
		}
	
	}

