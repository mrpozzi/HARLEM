harlem1Step<-function(fullData,Psi,deltaStar,tau,nGrid,lambda=Vectorize(function(x)1*((x>=10)&(x<=70))),Q0=NULL,T0=NULL,verbose=FALSE,display=FALSE){#,truth
	
	### Get the call and store the argument the function has been called with
	call <- match.call()
	argz <- lapply(2:length(call),function(i)eval(call[[i]]))
	names(argz) <- names(call)[2:length(names(call))]
	
	if(nGrid%%2==1)nGrid <- nGrid + 1
	xGrid <- seq(0,tau,length.out=nGrid); vGrid <- seq(0,tau,length.out=nGrid)
	
	w <- fullData$w
	
	obsData<-fullData[which(fullData$delta0*fullData$delta1==1),]
	#obsData <- obsData[which(obsData$x+deltaStar<=obsData$v),]
	obsData <- obsData[obsData$x!=0,]
	
	if(is.null(Q0)){
		if(verbose)cat("Initializing the Q's...\n")
		Q0 <- initQ(fullData,xGrid,vGrid,deltaStar=deltaStar,tau=tau,verbose=verbose)
		}
	Q1 <- Q0[[1]]; Q2 <- Q0[[2]]
	Sn <- Q1$Sn; Qc <- Q1$Qc; dGn <- Q1$dGn
	
	if(is.null(T0)){
		if(verbose)cat("Initializing the T's...\n")
		T0 <- initT(fullData,Qc,dGn, deltaStar=deltaStar,xGrid,vGrid)
		}
	T1 <- T0[[1]]; T2 <- T0[[2]]
	
	lambdaGrid <- lambda(xGrid); lambdaObs <- lambda(obsData$x)
	
	initPsi <- Psi(T1,T2,Q1,Q2,lambda,xGrid,vGrid)
	if(display|verbose)thetaHat <- Psi(T1,T2,Q1,Q2,lambda,xGrid,vGrid)
	
	
	obsData<-fullData[which(fullData$delta0*fullData$delta1==1),]
	obsData <- obsData[obsData$x!=0,]
	
	n <- nrow(obsData)
	
	slice1 <- lambdaGrid*apply(T1*Q2,2, .intFun,vGrid)
	slice2 <- lambdaGrid*apply(T2*Q2,2, .intFun,vGrid)
	grad <- t(apply(T1,1,function(t1){
		r <- t1/slice1
		r[(abs(t1)<=sqrt(.Machine$double.eps)&abs(slice1)<=sqrt(.Machine$double.eps))] <- 1
		r[lambdaGrid==0] <- 0
		return(r)
		})-apply(T2,1,function(t2){
			r <- t2/slice2
			r[(abs(t2)<=sqrt(.Machine$double.eps)&abs(slice2)<=sqrt(.Machine$double.eps))] <- 1
			r[lambdaGrid==0] <- 0
			return(r)
			}))/(n/nrow(fullData))
	grad[Q2==0] <- 0
	IC <- .linInt(grad,xGrid,vGrid)
	
	if(length(initPsi)==2){
		lambdax <- function(x) x*lambda(x);	lambdax2 <- function(x) x^2*lambda(x)
		i0lambda<-integrate(lambda,0,tau)$value
		i1lambda<-integrate(lambdax,0,tau)$value
		i2lambda<-integrate(lambdax2,0,tau)$value
		Lambda<-matrix(c(i0lambda,i1lambda,i1lambda,i2lambda),2,2)
		invLambda <- solve(Lambda)
		
		return(initPsi - invLambda %*% colSums(apply(cbind(1,obsData$x),2,"*",lambda(obsData$x)*IC(obsData$x, obsData$v)))/n)
		
		}else{
			invLambda <- 1/integrate(lambda,0,tau)$value
			return(initPsi - invLambda * sum(lambda(obsData$x)*IC(obsData$x, obsData$v))/n)
			}

	
	
	# return(new("HARLOT",thetaHat=cbind(initPsi,theta),
	                    # Q1=Q1,
	                    # Q2=Q2,
	                    # T1=T1,
	                    # T2=T2,
	                    # logL =logL,
	                    # normGrad=normGrad,
	                    # #envir="environment",
	                    # iter=1L,
	                    # call=call,
	                    # args=argz))
	
	}