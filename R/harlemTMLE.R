harlemTMLE<-function(fullData,targetStep,Psi, deltaStar,tau,nGrid,lambda=Vectorize(function(x)1*((x>=10)&(x<=70))),Q0=NULL,T0=NULL,verbose=FALSE,display=FALSE){#,truth
	
	
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
	
	if(display){
		dev.set(which=1)
		dev.set(which=1)
		dev.set(which=1)
		#dev.set(which=1)
		}
	
	count <- 0; cond <- FALSE
	targetStep <- targetStep(T1,T2,Q1,xGrid,vGrid,verbose,display)
	
	if(display|verbose)thetaHat <- Psi(T1,T2,Q1,Q2,lambda,xGrid,vGrid)
	
	epsOpt <- normGrad <- logL <- c()
	
	while(!cond){
		
		count<-count+1
		
		#if(verbose)
		cat(count," ")
		new <- targetStep(Q2)
		if(verbose)cat("\n")
		
		Q2 <- new$Q2; cond <- new$cond
		
		
		logL <- c(logL,new$logL)
		normGrad <-rbind(normGrad,new$normGrad)
		epsOpt <-rbind(epsOpt,new$epsOpt)
		
		#if(display|verbose) thetaHat<-rbind(thetaHat,Psi(T1,T2,Q1,Q2,lambda,intFun=intFun))
		#if(verbose)cat(" [",thetaHat[count,1],", ",thetaHat[count,2],"]\n")
		#if(verbose)cat("\n")
		
		if(display){
			
			dev.set(which=2)
			par(mfrow=c(1,2))
			plot(1:count, normGrad[,1],xlab="iter",ylab=expression(D[1]/sigma[1]),type="l",col="blue")
			abline(h=0,col="green")
			points(1:count, normGrad[,1],col="red",pch=16)
			plot(1:count, normGrad[,2],xlab="iter",ylab=expression(D[2]/sigma[2]),type="l",col="blue")
			points(1:count, normGrad[,2],col="red",pch=16)
			abline(h=0,col="green")
			# plot(0:count,thetaHat[,1],xlab="iter",ylab=expression(theta[1]),type="l",col="blue")
			# abline(h=truth[1],col="green")
			# points(0:count,thetaHat[,1],col=c("green",rep("red",count)),pch=16)
			# plot(0:count,thetaHat[,2],xlab="iter",ylab=expression(theta[2]),type="l",col="blue")
			# abline(h=truth[2],col="green")
			# points(0:count,thetaHat[,2],col=c("green",rep("red",count)),pch=16)	
			
			dev.set(which=3)
			plot(1:count,logL,xlab="iter",ylab="logLikelihood",type="l",col="blue")
			points(1:count,logL,col="red",pch=16)	
			
			}
			
			#if(count>1) cond <- abs(logLik[count]-logLik[count-1])<1e-5
		
		}
	
	
	cat("\n")
	theta <- Psi(T1,T2,Q1,Q2,lambda,xGrid,vGrid)
	
	return(list("theta"=theta,"Q"=list("Q1"=Q1,"Q2"=Q2),"T"=T0,"logL"=logL,"normGrad"= normGrad,"epsOpt"=epsOpt))
	
	}
