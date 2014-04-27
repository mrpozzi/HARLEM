initQ<-function(fullData,xGrid,vGrid, deltaStar,tau=100,tol=1e-3,k=10L,N=100L,verbose=FALSE,h=NULL){
	
	storage.mode(N) <- "integer"
	
	obsData <- fullData[which(fullData$delta0*fullData$delta1==1),]
	obsData <- obsData[obsData$x!=0,]
	# obsData <- obsData[which(obsData$x+deltaStar<=obsData$v),]
	
	QcFit<-survfit(Surv(fullData$v-fullData$w,1-fullData$delta1)~1)
	Qc<-stepfun(QcFit$time,c(1,QcFit$surv),f=0)
	Snfit<-survfit(Surv(fullData$w,fullData$v,fullData$delta1)~1)
	Sn<-stepfun(Snfit$time,c(1,Snfit$surv),f=0)
	dGn<-1/Sn(fullData$w); 
	dGn[Sn(fullData$w)==0] <- 0
	dGn[dGn/sum(dGn)>100/nrow(fullData)] <- 0
	dGn<-dGn/sum(dGn)
	
	if(Qc(deltaStar)<=0)stop("Initialization NOT well defined.\a")
	
	w <- fullData$w
	
	delta <- deltaStar
	nnn <- as.integer(nrow(obsData))
	vv <- obsData$v; xx <- obsData$x
	storage.mode(vv) <- storage.mode(xx) <- "double"
	nx <- length(xGrid); nv <- length(vGrid)
	Q2 <- matrix(1.0,nrow=nv,ncol=nx)#*1.0
	storage.mode(tau) <- "double"
	
	#step <- tau/N
	
	W <- c(1,2^((1:(N-1))%%2+1),1)
	W <- W%*%t(W)
	storage.mode(W) <- "double"
	
	ab <- range(log(xx),na.rm=TRUE); storage.mode(ab) <- "double"
	cd <- range(log(vv-xx-delta),na.rm=TRUE); storage.mode(cd) <- "double"
	
	cat("a=",ab[1]," b=",ab[2],"\n")
	cat("c=",cd[1]," d=",cd[2],"\n")
	
	if(is.null(h)){
		
		hOpt <- rep(0.0,2)
		
		h1Grid <- seq(0.05,25,length.out=k)
		h2Grid <- seq(0.05,25,length.out=k)
		delta0 <- diff(h1Grid)[1]
		
		# h1Grid<-seq(max(1.5,bandwidth.nrd(xx-vv)-5),bandwidth.nrd(xx-vv)+5,by=0.5)
		# h2Grid<-seq(max(1.5,bandwidth.nrd(vv)-3),bandwidth.nrd(vv)+5,by=0.5)
		
		n1<-length(h1Grid); n2<-length(h2Grid)
		
		storage.mode(h1Grid) <- "double"
		storage.mode(h2Grid) <- "double"
		
		m<-1
		if(verbose)cat("\n")
		while(!2^m*delta0/(k^m)<tol){
			deltam <- diff(h1Grid)[1]
			z <- .C('initQ', a_delta=delta,a_ab=ab, a_cd=cd,a_n=nnn,a_xx=xx,a_vv=vv,a_h1=h1Grid,a_h2=h2Grid,a_n1=n1,a_n2=n2,a_N=N,a_xGrid=xGrid,a_vGrid=vGrid,a_nx=nx,a_nv=nv,a_Q2=Q2, a_hOpt= hOpt, a_Weights=W)
			hOpt <- z$a_hOpt
			if(verbose)cat("*")
			
			h1Grid <- seq(max(0.1,hOpt[1]-deltam-tol), hOpt[1]+deltam+tol,length.out=k)
			h2Grid <- seq(max(0.1,hOpt[2]-deltam-tol), hOpt[2]+deltam+tol,length.out=k)
			
			# cat(m," ")
			# print(range(h1Grid))
			# print(range(h2Grid))
			
			m <- m+1
			
			}
		if(verbose)cat("\n")
		
		}else{
			
			hOpt <- rep(0.0,2)
			
			h1Grid <- h[1]
			h2Grid <- h[2]
			n1<-1L; n2<-1L
			
			storage.mode(h1Grid) <- "double"
			storage.mode(h2Grid) <- "double"
			
			z <- .C('initQ', a_delta=delta,a_ab=ab, a_cd=cd,a_n=nnn,a_xx=xx,a_vv=vv,a_h1=h1Grid,a_h2=h2Grid,a_n1=n1,a_n2=n2,a_N=N,a_xGrid=xGrid,a_vGrid=vGrid,a_nx=nx,a_nv=nv,a_Q2=Q2, a_hOpt=hOpt, a_Weights=W)
			
			
			}
		
		
		Q2<-z$a_Q2
		hOpt <- z$a_hOpt
	
	if(verbose)cat("h Opt.")
	print(hOpt)
	attr(Q2,"hOpt") <- hOpt
	#Q2[1,] <- Q2[,1] <- 0.0
	
	normC <- -1.0
	z <- .C('Norm', a_normC=normC,a_xGrid=xGrid,a_vGrid=vGrid,a_nx=nx,a_nv=nv,a_Fun=Q2)
	Q2 <- Q2/z$a_normC
	
	
	return(list("Q1"=list("Sn"=Sn,"dGn"=dGn,"Qc"=Qc),"Q2"=Q2))
	
	}

initT<-function(fullData,Qc,dGn, deltaStar, xGrid,vGrid){

	w <- fullData$w
	dGn <- dGn[order(w)]
	w <- sort(w)
	n <- length(w)
	storage.mode(w) <- "double"
	
	nx <- length(xGrid); nv <- length(vGrid)
	
	x <- get("x",envir=environment(Qc))
	y <- get("y",envir=environment(Qc))
	m <- length(x)
	
	T1 <- matrix(0.0,nrow=nv,ncol=nx)
	T2 <- matrix(0.0,nrow=nv,ncol=nx)
	
	z <- .C('initT', a_T1=T1, a_T2=T2, a_Stimes=x, a_Sjumps=y, a_dGn=dGn,a_delta= deltaStar,a_n=n,a_m=m,a_ww=w,a_xGrid=xGrid,a_vGrid=vGrid,a_nx=nx,a_nv=nv)
	
	T1 <- z$a_T1
	T2 <- z$a_T2
	
	
	return(list("T1"=z$a_T1,"T2"=z$a_T2))
	
	}
