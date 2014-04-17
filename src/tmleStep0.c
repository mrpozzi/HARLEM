#include <stdlib.h> 
#include <math.h>


#include "harlem.h"



void tmleStep0(double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv,double *a_lambda, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv,double *a_T1,double *a_T2, double *a_Q2,double *a_normGrad, double *a_logL, double *a_epsOpt){

  
  double slice1[*a_nx],slice2[*a_nx];
  double hv=(a_vGrid[*a_nv-1]-a_vGrid[0])/2.0;
  ///printf("%g\n",hx);
  //printf("\n");

  /// Compute the gradient on the grid
  for(int j=0;j<*a_nx;j++)
    {

      slice1[j]=a_T1[*a_nv*j]*a_Q2[*a_nv*j]+a_T1[*a_nv-1+*a_nv*j]*a_Q2[*a_nv-1+*a_nv*j];
      slice2[j]=a_T2[*a_nv*j]*a_Q2[*a_nv*j]+a_T2[*a_nv-1+*a_nv*j]*a_Q2[*a_nv-1+*a_nv*j];

      for(int k=1;k<*a_nv/2;k++)
	{
	  slice1[j]+=2*a_T1[2*k+*a_nv*j]*a_Q2[2*k+*a_nv*j]+4*a_T1[2*k-1+*a_nv*j]*a_Q2[2*k-1+*a_nv*j];
	  slice2[j]+=2*a_T2[2*k+*a_nv*j]*a_Q2[2*k+*a_nv*j]+4*a_T2[2*k-1+*a_nv*j]*a_Q2[2*k-1+*a_nv*j];
	}

      slice1[j]*=hv/3.0;
      slice2[j]*=hv/3.0;

    }

  /// Compute the Gradient and the boundaries for the optimization problem
  double grad[*a_nx*(*a_nv)];
  double gradMax=0.0;

  for(int i=0;i<*a_nv;i++)
    {
      for(int j=0;j<*a_nx;j++)
	{

	    if(ABS(slice1[j])<sqrt(DBL_EPSILON)||ABS(slice2[j])<sqrt(DBL_EPSILON))
	    {
	      grad[i+*a_nv*j]=0.0;
	    }else
	    {
	      grad[i+*a_nv*j]=(a_T1[i+*a_nv*j]/slice1[j]-a_T2[i+*a_nv*j]/slice2[j]);
	    }


	  grad[i+*a_nv*j]= a_lambda[j]*grad[i+*a_nv*j];
	  
	  double absGrad=ABS(grad[i+*a_nv*j]);

	  if(absGrad>gradMax) gradMax=absGrad;

	}
    }
  double limEpsx = 1/(gradMax);

  /// Evaluate the Gradient @ the observed data
  double gradObs[*a_n];
  Interp(a_xx,a_vv,a_n,gradObs,a_xGrid,a_vGrid,a_nx,a_nv,grad);

  /// Solve the Optimization Problem and find the epsilon
  double epsOpt = 0.0;
  double fMax = critFun0(a_xx,gradObs,*a_n,epsOpt,gradMax);

  double step1=2.0*limEpsx/(nEps-1);

  //printf("step=[%g,%g]\n",step1,step2);

  double eps=-limEpsx;
  for(int i=0;i<nEps;i++)
  {
    
    double f = critFun0(a_xx,gradObs,*a_n,eps,gradMax);
	 
    if(f>fMax)
      {
	fMax=f;
	epsOpt=eps;
      }
    eps=MIN(eps+step1,limEpsx);
	
  }

  *a_epsOpt=epsOpt;

  /// Compute Q2 @ observed data + compute empirical Mean and Variance
  double empMean=0.0;
  double empVar=0.0;
  double sumsq=0.0;
  double offSet = 0.0;

  for(int i=0;i<*a_n;i++)
    {
      int ii = 1;
      double q2Obs;
      Interp(&a_xx[i],&a_vv[i],&ii,&q2Obs,a_xGrid,a_vGrid,a_nx,a_nv,a_Q2);
      offSet += log(q2Obs);
 
      empMean += (gradObs[i] - gradObs[0]);
      sumsq += (gradObs[i] - gradObs[0])*(gradObs[i] - gradObs[0]);

    }
  offSet/=*a_n;
  *a_logL=offSet + fMax;

  empVar = ((sumsq - empMean*empMean/(*a_n))/(*a_n-1))/(*a_n);
  empMean /= *a_n;
  empMean += gradObs[0];

  *a_normGrad=empMean/sqrt(empVar);
  
  /// Update Q2
  for(int i=0;i<*a_nv;i++)
    {
      for(int j=0;j<*a_nx;j++)
	{
	  a_Q2[i+*a_nv*j]=(1+grad[i+*a_nv*j]*epsOpt)*a_Q2[i+*a_nv*j];
	  ///if(a_Q2[i+*a_nv*j]<0.0)printf("(%g,%g) %g %g (%g)\n",eps1Opt,eps2Opt,grad[i+*a_nv*j],(eps1Opt+eps2Opt*a_xGrid[j]),a_xGrid[j]);
	}
    }
  
  /// Normalize Q2
  double normC;
  Norm(&normC,a_xGrid,a_vGrid,a_nx,a_nv, a_Q2);

  //printf("\n%g\n",normC);

  for(int i=0;i<*a_nv;i++)
    {
      for(int j=0;j<*a_nx;j++)
	{
	  a_Q2[i+*a_nv*j]/=normC;
	}
    }

}




double critFun0(double *a_xx, double *a_gradObs,int a_n,double a_eps, double a_D){

  double val = 0.0;
  if(ABS(a_eps)*a_D>1)return(-1.0);
  for(int i=0;i<a_n;i++)
    {      
      val += log(1.0+a_eps*a_gradObs[i]);
    }
  val/=a_n;
  return(val);
  
}
