.intFun <- function(f,nodes){
	
	if(length(f)!=length(nodes))error("Number of nodes different from number of function Values!!!")
	
	ab <- range(nodes); m <- length(nodes)
	h <- (ab[2]-ab[1])/m
	
	return(h/3*(f[1]+2*sum(f[(1:length(f))[1:length(f)%%2==0]])+4*sum(f[(1:length(f))[1:length(f)%%2==1]])+f[length(f)]))
	
	}

	
	
.linInt <- function(fVals,xVals,yVals){
	
	nx<-length(xVals); ny<-length(yVals)
	
	return(function(x,y){
		

		n<-length(x)
		phi <- rep(0.0,n)	

		z <- .C('Interp', a_x=x, a_y=y, a_n=n, a_Val=phi,a_xGrid= xVals,a_vGrid= yVals,a_nx=nx,a_nv=ny,a_Fun= fVals)
		return(z$a_Val)
		
		})
	
	}