#include <stdlib.h> 
#include <math.h>

#include "harlem.h"



/// Computes the normalizing constant of a function by linear interpolation
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

	  *a_normC+=(a_xGrid[j+1]-a_xGrid[j])*(a_vGrid[i+1]-a_vGrid[i])*((phi00+phi01+phi11)/3+(phi00+phi10+phi11)/3)/2;

	}
    }
}


/// linear interpolation
void Interp(double *a_x,double *a_y, int *a_n, double *a_Val, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun){


  for(int k=0;k<*a_n;k++)
    {

      int i=0;int i0;
      int found=0;

      while(found==0&&i<*a_nv-1)
	{
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
	  if(a_xGrid[j]>a_x[k])
	    {
	      found=1;
	      j0=j-1;
	    }
	  j++;
	}
      
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
  



int Cholesky(double* A, long n){
  
  int i,j,k;
  double sum;
  
  /*Cholesky Decomposition!!!*/
  for (i=0; i<n; i++) {
	
	sum = 0;
	for (k=0; k<=i-1; k++)
	  sum += A[n*i+k]*A[n*i+k];
	
	if(*(A+(n+1)*i)-sum <= 0){//check for rounding errors
	    //printf("\nERROR: Cholesky Failed!!!\a\n");
	    exit(1);
	}else{
	    A[(n+1)*i] = sqrt(A[(n+1)*i]-sum);
	}
	
	for (j=i+1; j<n; j++) {
	    for (k=i-1,sum=*(A+n*j+i); k>=0; k--)
		sum -= *(A+n*i+k) * *(A+n*j+k);
	    *(A+n*j+i) = sum/(*(A+(n+1)*i));
	}	
  }
  
  for (i=0; i<n; i++) 
      for (j=0; j<i; j++) 
	  *(A+n*j+i)=0; 
  
  return(0);
  
  
}


int UGauss(double* A,double* b,long n) {
  
  int i, j;
  double sum;
  
  *(b+n-1) /= *(A+(n+1)*(n-1));
  for(i=n-2;i>=0;i--){
	sum = *(b+i);
	for(j=i+1;j<n;j++){
	  sum -= *(A+i+n*j) * *(b+j);
	}
	*(b+i) = sum / *(A+(n+1)*i);
  }
  
  return(0);
  
}

/// This routine takes an upper triangular matrix but treats it as a lower triangula one...
int LGauss(double* A,double* b,long n) {
  
  int i, j;
  double sum;
  
  for(i=0;i<n;i++){
	sum = *(b+i);
	for(j=0;j<i;j++){
	  sum -= *(A+j+n*i) * *(b+j);
	}
	*(b+i) = sum / *(A+(n+1)*i);
  }
  
  return(0);
  
}
