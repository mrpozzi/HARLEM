#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>

#include "harlem.h"



/// TMLE step for linear 
void tmleStep1(double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv,double *a_lambda, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv,double *a_T1,double *a_T2, double *a_Q2,double *a_normGrad, double *a_logL, double *a_epsOpt){

  
  double slice1[*a_nx],slice2[*a_nx];
  double hv=(a_vGrid[*a_nv-1]-a_vGrid[0])/(*a_nv);

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
  double xGradMax=0.0;

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


	  grad[i+*a_nv*j] *= a_lambda[j];
	  
	  double absGrad = ABS(grad[i+*a_nv*j]);
	  double absXgrad = a_xGrid[j]*absGrad;

	  if(absGrad>gradMax) gradMax = absGrad;
	  if(absXgrad>xGradMax) xGradMax = absXgrad;

	}
    }

  double limEpsx = 1/(gradMax); double limEpsy = 1/(xGradMax);

  /// Evaluate the Gradient @ the observed data
  double gradObs[*a_n];
  Interp(a_xx,a_vv,a_n,gradObs,a_xGrid,a_vGrid,a_nx,a_nv,grad);


  /// Solve the Optimization Problem and find the epsilon
  double eps1Opt = 0.0;
  double eps2Opt = 0.0;
  double fMax = critFun1(a_xx,gradObs,*a_n,eps1Opt,eps2Opt,gradMax,xGradMax);

  double step1 = 2.0*limEpsx/(nEps-1);
  double step2 = 2.0*limEpsy/(nEps-1);


  double eps1 = -limEpsx;
  for(int i=0;i<nEps;i++)
  {
    
    double eps2 = -limEpsy;
    
    for(int j=0;j<nEps;j++)
      {
	
	double f = critFun1(a_xx,gradObs,*a_n,eps1,eps2,gradMax,xGradMax);
	 
	  if(f>fMax)
	    {
	      fMax = f;
	      eps1Opt = eps1;
	      eps2Opt = eps2;
	    }
	  eps2 = MIN(eps2+step2,limEpsy);

      }
    eps1 = MIN(eps1+step1,limEpsx);
	
  }

  a_epsOpt[0]=eps1Opt;
  a_epsOpt[1]=eps2Opt;


  /// Compute Q2 @ observed data + compute empirical Mean and Variance
  double empMean[2] = {0.0,0.0};
  double empVar[2] = {0.0,0.0};
  double sumsq[2] = {0.0,0.0};
  double offSet = 0.0;

  for(int i=0;i<*a_n;i++)
    {
      int ii = 1;
      double q2Obs;
      Interp(&a_xx[i],&a_vv[i],&ii,&q2Obs,a_xGrid,a_vGrid,a_nx,a_nv,a_Q2);
      offSet += log(q2Obs);
 
      empMean[0] += (gradObs[i] - gradObs[0]);
      empMean[1] += (a_xx[i]*gradObs[i] - a_xx[0]*gradObs[0]);

      sumsq[0] += (gradObs[i] - gradObs[0])*(gradObs[i] - gradObs[0]);
      sumsq[1] += (a_xx[i]*gradObs[i] - a_xx[0]*gradObs[0])*(a_xx[i]*gradObs[i] - a_xx[0]*gradObs[0]);

    }
  offSet/=*a_n;
  *a_logL=offSet + fMax;

  empVar[0] = ((sumsq[0] - empMean[0]*empMean[0]/(*a_n))/(*a_n-1))/(*a_n);
  empVar[1] = ((sumsq[1] - empMean[1]*empMean[1]/(*a_n))/(*a_n-1))/(*a_n);


  empMean[0] /= *a_n;
  empMean[1] /= *a_n;

  empMean[0] += gradObs[0];
  empMean[1] += a_xx[0]*gradObs[0];

  a_normGrad[0]=empMean[0]/sqrt(empVar[0]);
  a_normGrad[1]=empMean[1]/sqrt(empVar[1]);
  
  /// Update Q2
  for(int i=0;i<*a_nv;i++)
    {
      for(int j=0;j<*a_nx;j++)
	{
	  a_Q2[i+*a_nv*j]=(1+grad[i+*a_nv*j]*(eps1Opt+eps2Opt*a_xGrid[j]))*a_Q2[i+*a_nv*j];
	}
    }
  
  /// Normalize Q2
  double normC;
  Norm(&normC,a_xGrid,a_vGrid,a_nx,a_nv, a_Q2);

  for(int i=0;i<*a_nv;i++)
    {
      for(int j=0;j<*a_nx;j++)
	{
	  a_Q2[i+*a_nv*j]/=normC;
	}
    }

}




double critFun1(double *a_xx, double *a_gradObs,int a_n,double a_eps1, double a_eps2, double a_D, double a_xD){

    if(ABS(a_eps1)*a_D+ABS(a_eps2)*a_xD>1)return(-1.0);
    double val = 0.0;
    for(int i=0;i<a_n;i++)
    {
	val += log(1.0+(a_eps1+a_eps2*a_xx[i])*a_gradObs[i]);
    }
    val/=a_n;
    return(val);
  
}
