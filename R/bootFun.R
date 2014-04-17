bootHARLEM <- function(harlem,B=1000,deltaFun=function(data)max(diff(sort(data$w)))){
	
	fullData<-harlem@args$fullData
	lambda <- harlem@args$lambda
	nGrid <- harlem@args$nGrid
	tau <- harlem@args$tau
	Q1 <- harlem@Q1
	Q2tmle <- harlem@Q2
	T1 <- harlem@T1
	T2 <- harlem@T2
	targetStep <- harlem@args$targetStep
	Psi <- harlem@args$Psi

	xGrid <- seq(0,tau,length.out=nGrid); vGrid <- seq(0,tau,length.out=nGrid)
	
	tmlePsi <- harlem@thetaHat[,2]
	
	if(length(tmlePsi)==2){
		
		lambdax <- function(x) x*lambda(x);	lambdax2 <- function(x) x^2*lambda(x)
		i0lambda<-integrate(lambda,0,tau)$value
		i1lambda<-integrate(lambdax,0,tau)$value
		i2lambda<-integrate(lambdax2,0,tau)$value
		Lambda<-matrix(c(i0lambda,i1lambda,i1lambda,i2lambda),2,2)
		invLambda <- solve(Lambda)
		}else{
			invLambda <- 1/integrate(lambda,0,tau)$value
			}
	
	lambdaGrid <- lambda(xGrid)
	
	obsData<-fullData[which(fullData$delta0*fullData$delta1==1),]
	obsData <- obsData[obsData$x!=0,]
	deltaStar <- deltaFun(obsData)
	
	n <- nrow(obsData)
	
	slice1 <- lambdaGrid*apply(T1*Q2tmle,2, .intFun,vGrid)
	slice2 <- lambdaGrid*apply(T2*Q2tmle,2, .intFun,vGrid)
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
			}))
	IC <- .linInt(grad,xGrid,vGrid)
	
	bootFun <- function(b){
		
		cat(paste(100*b/B,"% ",sep=""))
		
		bootData <- fullData[sample(1:nrow(fullData),replace=TRUE),]
		obsData<-bootData[which(bootData$delta0*bootData$delta1==1),]
		
		lambdaObs <- lambda(obsData$x)
		deltaStar <- deltaFun(obsData)
		
		QcFit<-survfit(Surv(bootData$v-bootData$w,1-bootData$delta1)~1)
		Qc<-stepfun(QcFit$time,c(1,QcFit$surv),f=0)
		Snfit<-survfit(Surv(bootData$w,bootData$v,bootData$delta1)~1)
		Sn<-stepfun(Snfit$time,c(1,Snfit$surv),f=0)
		dGn<-1/Sn(bootData$w)
		dGn[Sn(bootData$w)==0] <- 0; dGn<-dGn/sum(dGn)
		bootQ1 <- list("Sn"=Sn,"dGn"=dGn,"Qc"=Qc)
		
		bootT0 <- initT(bootData,Qc,dGn,deltaStar,xGrid,vGrid); bootT1 <- bootT0[[1]]; bootT2 <- bootT0[[2]]
		
		assign("v",bootData$v,envir=environment(Psi))
		assign("w",bootData$w,envir=environment(Psi))
		assign("obsData",obsData,envir=environment(Psi))
		assign("fullData",bootData,envir=environment(Psi))
		bootPsi <- Psi(bootT1,bootT2,bootQ1,Q2tmle,lambda,xGrid,vGrid)
		
		if(length(tmlePsi)==2){
			return(bootPsi - tmlePsi - invLambda %*% colSums(apply(cbind(1,obsData$x),2,"*",lambda(obsData$x)*IC(obsData$x, obsData$v)))/n)
			}else{
				return(bootPsi - tmlePsi - invLambda * sum(lambda(obsData$x)*IC(obsData$x, obsData$v))/n)
				}
		
		
		}
			
	bootRep <- lapply(1:B,bootFun)
	
	cat("\n")
	
	#return(apply(do.call(cbind,bootRep),1,var))
	return(do.call(cbind,bootRep))
	
	}
