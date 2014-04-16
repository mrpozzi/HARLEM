#include <stdlib.h> 
#include <stdio.h> 
#include<math.h>



void Norm(double *a_normC, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun);
void Interp(double *a_x,double *a_y, int *a_n, double *a_Val, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun);

void Norm(double *a_normC, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun){

  *a_normC=0.0;

  for(int i=0;i<*a_nv-1;i++)
    {
      for(int j=0;j<*a_nx-1;j++)
	{
	  
	  double phi00=a_Fun[i+*a_nv*j];
	  double phi10=a_Fun[i+*a_nv*(j+1)];
	  double phi01=a_Fun[i+1+*a_nv*j];
	  double phi11=a_Fun[i+1+*a_nv*(j+1)];
	  
	  /*
	  double psi00=a_Fun[i+1+*a_nv*j];
	  double psi01=a_Fun[i+1+*a_nv*(j+1)];
	  double psi10=a_Fun[i+*a_nv*j];
	  double psi11=a_Fun[i+1+*a_nv*(j+1)];
	  */


	  *a_normC+=(a_xGrid[j+1]-a_xGrid[j])*(a_vGrid[i+1]-a_vGrid[i])*((phi00+phi01+phi11)/3+(phi00+phi10+phi11)/3)/2;

	}
    }
}



void Interp(double *a_x,double *a_y, int *a_n, double *a_Val, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun){


  for(int k=0;k<*a_n;k++)
    {

      int i=0;int i0;
      int found=0;


      /*
	int up=a_n-1;
	int low=0;
	int i = ceil((up+low)/2);
	while(found==0&&i<*a_nv-1)
	{
	
	if(a_Stimes[i]>a_val)
	{
	if(a_Stimes[i-1]<=a_val)
	    {
	      return(a_Sjumps[i-1]);
	    }else 
	    {
	      up=i;
	      i = ceil((up+low)/2);
	      }
	      }else if(a_Stimes[i+1]>=a_val)
	      {
	      return(a_Sjumps[i]);
	      }else 
	      {
	      
	      low=i;
	      i = ceil((up+low)/2);
	      }
	      }*/
      while(found==0&&i<*a_nv-1)
	{
	  //if(a_vGrid[i]<=a_y[k]&&a_vGrid[i+1]>a_y[k])
	  if(a_vGrid[i]>a_y[k])
	    {
	      found=1;
	      i0=i-1;
	    }
	  i++;
	}

      found=0;
      int j=0;
      int j0;
      while(found==0&&j<*a_nx-1)
	{
	  	  
	  //if(a_xGrid[j]<=a_x[k]&&a_xGrid[j+1]>a_x[k])
	  if(a_xGrid[j]>a_x[k])
	    {
	      found=1;
	      j0=j-1;
	    }
	  j++;
	}

      //printf("%d %d\n",i0,j0);
      /*
	# ix0<-max(which(xVals<=x)); iy0<-max(which(yVals<=y))
	# x0<-xVals[ix0]; x1<-xVals[ix0+1]; y0<-yVals[iy0]; y1<-yVals[iy0+1]
	# phi00<-fVals[ix0,iy0]; phi01<-fVals[ix0,iy0+1]
	# phi10<-fVals[ix0+1,iy0]; phi11<-fVals[ix0+1,iy0+1]
	# d <- as.numeric((y-y0)*(x1-x0)>=(x-x0)*(y1-y0))
	# phi <- phi00+(x-x0)/(x1-x0)*(d*(phi11-phi01)+(1-d)*(phi10-phi00))+(y-y0)/(y1-y0)*(d*(phi01-phi00)+(1-d)*(phi11-phi10))
      */
      
      /*
      double psi00=a_Fun[i0+1+*a_nv*j0];
      double psi01=a_Fun[i0+1+*a_nv*(j0+1)];
      double psi10=a_Fun[i0+*a_nv*j0];
      double psi11=a_Fun[i0+1+*a_nv*(j0+1)];
      */
      
      double phi00=a_Fun[i0+*a_nv*j0];
      double phi10=a_Fun[i0+*a_nv*(j0+1)];
      double phi01=a_Fun[i0+1+*a_nv*j0];
      double phi11=a_Fun[i0+1+*a_nv*(j0+1)];
      

      a_Val[k]=phi00;
      
      if((a_y[k]-a_vGrid[i0])*(a_xGrid[j0+1]-a_xGrid[j0])>=(a_vGrid[i0+1]-a_vGrid[i0])*(a_x[k]-a_xGrid[j0]))
	{
	  a_Val[k]+=((a_x[k]-a_xGrid[j0])/(a_xGrid[j0+1]-a_xGrid[j0]))*(phi11-phi01)+((a_y[k]-a_vGrid[i0])/(a_vGrid[i0+1]-a_vGrid[i0]))*(phi01-phi00);
	}else
	{
	  a_Val[k]+=((a_x[k]-a_xGrid[j0])/(a_xGrid[j0+1]-a_xGrid[j0]))*(phi10-phi00)+((a_y[k]-a_vGrid[i0])/(a_vGrid[i0+1]-a_vGrid[i0]))*(phi11-phi10);
	}
      
    }
}
