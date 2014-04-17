harlemTMLE<-function(fullData,targetStep,Psi, deltaStar,tau,nGrid,lambda=Vectorize(function(x)1*((x>=10)&(x<=70))),initPsi=NULL,Q0=NULL,T0=NULL,verbose=FALSE,display=FALSE){#,truth
	
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
	
	if(display){
		dev.set(which=1)
		dev.set(which=1)
		dev.set(which=1)
		#dev.set(which=1)
		}
	
	count <- 0L; cond <- FALSE
	targetStep <- targetStep(T1,T2,Q1,xGrid,vGrid,verbose,display)
	# argz[["targetStep"]] <- targetStep
	
	initPsi <- Psi(T1,T2,Q1,Q2,lambda,xGrid,vGrid)
	if(display|verbose)thetaHat <- Psi(T1,T2,Q1,Q2,lambda,xGrid,vGrid)
	
	normGrad <- logL <- c()
	
	while(!cond){
		
		count<-count+1L
		
		#if(verbose)
		cat(count," ")
		new <- targetStep(Q2)
		if(verbose)cat("\n")
		
		Q2 <- new$Q2; cond <- new$cond
		
		
		#epsOpt <-rbind(epsOpt,new$epsOpt)
		
		#if(display|verbose) thetaHat<-rbind(thetaHat,Psi(T1,T2,Q1,Q2,lambda,intFun=intFun))
		#if(verbose)cat(" [",thetaHat[count,1],", ",thetaHat[count,2],"]\n")
		#if(verbose)cat("\n")
		
		if(display){
			
			logL <- c(logL,new$logL)
			normGrad <-rbind(normGrad,new$normGrad)
		
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
			
			}else{
				logL <- new$logL
				normGrad <- as.matrix(new$normGrad)
				}
			
			#if(count>1) cond <- abs(logLik[count]-logLik[count-1])<1e-5
		
		}
	
	
	cat("\n")
	theta <- Psi(T1,T2,Q1,Q2,lambda,xGrid,vGrid)
	
	# return(list("theta"=theta,"Q"=list("Q1"=Q1,"Q2"=Q2),"T"=T0,"logL"=logL,"normGrad"= normGrad,"epsOpt"=epsOpt))
	
	return(new("HARLOT",thetaHat=cbind(initPsi,theta),
	                    Q1=Q1,
	                    Q2=Q2,
	                    T1=T1,
	                    T2=T2,
	                    logL =logL,
	                    normGrad=normGrad,
	                    #envir="environment",
	                    iter=count,
	                    call=call,
	                    args=argz))
	
	}


#################################
### H Accelerated Residual Lifetime ObjecT

setClass("HARLOT",representation=representation(thetaHat="matrix",
                                             Q1="list",
                                             Q2="matrix",
                                             T1="matrix",
                                             T2="matrix",
                                             logL = "numeric",
                                             normGrad ="matrix",
                                             #envir="environment",
                                             iter="integer",
                                             call="call",
                                             args="list"))


### Minimal description
setMethod("show",signature("HARLOT"),function(object){
	cat("Call:\n")
	print(object@call)
	#type <- 
	cat("Psi:\n")
	print(object@thetaHat[,2])
	})

### more extensive summary...
setMethod("summary",signature("HARLOT"),function(object){
	cat("Call:\n")
	print(object@call)
	cat("Psi TMLE:\n")
	print(object@thetaHat[,2])
	cat("(",object@iter,"Iterations)\n")
	cat("Log Likelihood: ")
	cat(object@logL[length(object@logL)],"\n")
	cat("Gradient: ")
	print(object@normGrad[nrow(object@normGrad),])
	cat("Initial Psi:\n")
	print(object@thetaHat[,1])
	
	})

### plot method
setMethod("plot", signature = "HARLOT", definition = function (x,whichOne=1L:3L,xlim=c(-30,30),ylim=c(-30,30),main,phi=30, theta=-30,border=NA,xlab=expression(x[1]),ylab=expression(x[2]),zlab="target",...)  {
	
	show <- rep(FALSE, 5)
    show[whichOne] <- TRUE
	
	if(length(whichOne)>1) par(ask=TRUE)
	if(missing(main)) main <- paste(x@call[[1]],"Algorithm;",x@call[["target"]],"Target",sep=" ")
	
	Q2 <- x@Q2
	T1<-x@T1;T2<-x@T2
	lambda <- x@args$lambda
	nGrid<-x@args$nGrid
	if(nGrid%%2==1)nGrid <- nGrid + 1
	xGrid <- seq(0,x@args$tau,length.out=nGrid)
	vGrid <- seq(0,x@args$tau,length.out=nGrid)
	lambdaGrid <- lambda(xGrid)
	
			
	###	3-D plot
	if(show[1L]) {
		
		v <- x@args$fullData$v	
		
		
		Sn <- x@Q1$Sn
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
		
		plot(xGrid[lambdaGrid!=0],R[lambdaGrid!=0],type="l",col="green",ylab="Psi",xlab="x")
		
		
		linFun <- function(xx){x@thetaHat[1,1]+x@thetaHat[2,1]*xx}
		linFun<-Vectorize(linFun)
		plot(linFun,xlim=c(65,105),col="red",add=TRUE)
		linFun <- function(xx){x@thetaHat[1,2]+x@thetaHat[2,2]*xx}
		linFun<-Vectorize(linFun)
		plot(linFun,xlim=c(65,105),col="blue",add=TRUE)
		
		legend("bottomleft",legend=c(expression(Psi[N]),expression(psi[TMLE]),expression(psi[N])),col=c("green","blue","red"),lty=1)


    	} 
    
    ### contour
    if(show[2L]&x@args$display){    
    	normGrad <- x@normGrad
    	count<-x@iter
    	par(mfrow=c(1,2))
    	plot(1:count, normGrad[,1],xlab="iter",ylab=expression(D[1]/sigma[1]),type="l",col="blue")
    	abline(h=0,col="green")
    	points(1:count, normGrad[,1],col="red",pch=16)
    	plot(1:count, normGrad[,2],xlab="iter",ylab=expression(D[2]/sigma[2]),type="l",col="blue")
    	points(1:count, normGrad[,2],col="red",pch=16)
    	abline(h=0,col="green")
    	

    	}
    		   
    if(show[3L]&x@args$display){
    	logL <- x@logL
    	count<-x@iter
    	par(mfrow=c(1,1))
		plot(1:count,logL,xlab="iter",ylab="logLikelihood",type="l",col="blue")
			points(1:count,logL,col="red",pch=16)
		}
		
	if(show[3L]){
		w <- x@args$fullData$w
		Sn <- x@Q1$Sn
		Qc <- x@Q1$Qc
		dGn <- x@Q1$dGn
		par(mfrow=c(1,3))
		plot(Sn,main=expression(S[n]))
		plot(Qc,main=expression(Q[c]))
		Gn <- function(ww){sum(dGn[order(w[w<ww])])}
		Gn <- Vectorize(Gn)
		plot(Gn,main=expression(G[n]),xlim=range(w))
		}
							
	
	})
	
### subsetting method
setMethod("[[", signature = "HARLOT", definition = function (x, i,j,..., drop = TRUE) {
	if (!missing(j)) {stop("Wrong number of dimensions.")}
	if (!missing(i)) {
		return(switch(class(i),"character" = attributes(x)[[i]],
		                       "integer" = attributes(x)[[i]],
		                       "numeric" = attributes(x)[[i]],
		                       "logical" = attributes(x)[c(i,rep(FALSE,length(attributes(x))-length(i)))],
		                        stop("Subsetting object needs to be either a character, a numeric/integer or a logical.")
		                        ))
		}else{return(NULL)}
  
   })

### subsetting method
setMethod("$", signature = "HARLOT", definition = function(x, name) {
	 x[[name]]
	 })

### alias for slotNames 
setMethod("names", signature = "HARLOT", definition = function(x) {
	 slotNames(x)
	 })
	 

# setMethod("attach",signature("harlemTMLE"),function(what, warn.conflicts = TRUE){#,...
	# attach(what@envir, warn.conflicts= warn.conflicts)#,...
	# })

# setMethod("detach",signature("harlemTMLE"),function(name){
	# detach(name@envir)
	# })
