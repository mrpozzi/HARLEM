stepConst <- function(fullData,deltaStar,tau,lambda,tol){

	obsData <- fullData[which(fullData$delta0*fullData$delta1==1),]
	obsData <- obsData[obsData$x!=0,]
	#obsData <- obsData[which(obsData$x+deltaStar<=obsData$v),]
	
	w <- fullData$w; v <- obsData$v
	x <- obsData$x; n <- nrow(obsData)
	
	function(T1,T2,Q1,xGrid,vGrid,verbose=FALSE,display=FALSE){
		
		lambdaGrid <- lambda(xGrid)
		nx <- length(xGrid); nv <- length(vGrid)
		
		function(Q2){
			
			delta<- deltaStar
			storage.mode(lambdaGrid)<-"double"
			
			logL <- 0.0
			normGrad <- 0.0
			epsOpt <- 0.0
			
			z <- .C('tmleStep0',a_delta=delta,a_tau=tau,a_n=n,a_xx=x,a_vv=v,a_lambda=lambdaGrid, a_xGrid=xGrid,a_vGrid=vGrid,a_nx=nx,a_nv=nv,a_T1=T1,a_T2=T2,a_Q2=Q2, a_normGrad=normGrad,a_logL=logL,a_epsOpt=epsOpt)
			
			Q2 <- z$a_Q2
			normGrad <- z$a_normGrad
			logL <- z$a_logL
			epsOpt <- z$a_epsOpt
			
			if(any(Q2<0))stop("Negative Elements in Q2...\a")
			
			
			if(verbose) cat("D: ", normGrad," |  logLik =",logL," ")
			
			cond <- (max(abs(normGrad))<tol)
			
			return(list("Q2"=Q2,"cond"=cond,"normGrad"= normGrad,"logL"=logL,"epsOpt"=epsOpt))
			
			}
		
		}
	}

constPsi <- function(fullData,deltaStar,tau){
	
	w <- fullData$w; v <- fullData$v
	
	obsData <- fullData[which(fullData$delta0*fullData$delta1==1),]
	obsData <- obsData[obsData$x!=0,]
	# obsData <- obsData[which(obsData$x+deltaStar<=obsData$v),]
	
	function(T1,T2,Q1,Q2,lambda,xGrid,vGrid){
		
		Sn <- Q1$Sn; Qc <- Q1$Qc; dGn <- Q1$dGn
		
		
		
		Lambda<-integrate(lambda,0,tau)$value
		
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

		intR <- .intFun(R*lambdaGrid,xGrid)
		
		theta <- intR/Lambda

		return(theta)
		
		}
	
	}

