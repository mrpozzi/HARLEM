#include <stdio.h>

#include <stdlib.h>
#include <math.h>

#include "harlem.h"



/// TMLE step 
void tmleStep(double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv,double *a_lambda, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv,double *a_T1,double *a_T2, double *a_Q2,double *a_normGrad, double *a_logL, double *a_epsOpt, int *a_degree){

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

  /// Compute the Gradient 
  double grad[*a_nx*(*a_nv)];

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

	}
    }

  /// Evaluate the Gradient @ the observed data
  double gradObs[*a_n];
  Interp(a_xx,a_vv,a_n,gradObs,a_xGrid,a_vGrid,a_nx,a_nv,grad);


  /// Solve the Optimization Problem and find the epsilon
  double gradf[*a_degree+1];
  double delta[*a_degree+1];
  double Hess[(*a_degree+1)*(*a_degree+1)];

  double decrement = 10*bigN;
  double decrementOld=0.0;
  double eps[*a_degree+1];  
      
  
//  double fOld=0.0;


  //while (MIN(decrement/2,ABS(decrement-decrementOld)) >= sqrt(DBL_EPSILON))
  double m = (double)*a_nv * (double)*a_nv;
  double barr = bigN;

  printf("m = %lf : barr = %lf\n",m,barr);

  double fNew = critFun(a_epsOpt, a_xx, gradObs, *a_n, *a_degree, gradf, Hess, grad, a_xGrid, *a_nv, *a_nx,barr);

  int Iter = 0;

  while(1){

      Iter++;
      int iter = 0;
 
      while (1)
      {
	  iter++;
	  printf("[%d](%d) ",Iter,iter);

	  printf("f(theta)=%g | ",fNew);

	  /* solve the system */
	  Cholesky(Hess,*a_degree+1);
	  
	  LGauss(Hess, gradf, *a_degree+1); ///1
	  
	  decrement = 0.0;
	  for(int k=0;k<=*a_degree;k++)
	  {
	      decrement += gradf[k]*gradf[k];
	  }
	  
	  UGauss(Hess, gradf, *a_degree+1); ///2

	  printf("%lf \n",decrement);
	  if(decrement/2 < sqrt(DBL_EPSILON) && iter > 1) break;

	  double fOld = fNew;
	  double t = 1;

	  for(int k=0;k<=*a_degree;k++) delta[k] = gradf[k];

	  do{
 
	      for(int k=0;k<=*a_degree;k++) eps[k] = a_epsOpt[k] - t*delta[k];

	      fNew = critFun(eps, a_xx, gradObs, *a_n, *a_degree, gradf, Hess, grad, a_xGrid, *a_nv, *a_nx, barr);
	      t *= beta;
	      //printf("*");
	      
	  } while (fNew > fOld-alpha*t*decrement);
      
	  for(int k=0;k<=*a_degree;k++) a_epsOpt[k] = eps[k];
	  printf("[t=%g | f=%g]\n",t/beta,fNew);
      
      }
      printf("%lf\n",m/barr);
      if(m/barr < sqrt(DBL_EPSILON)) break;
      barr *= mu;
      //break;
  } 

  printf("\n");
  printf("eps: %lf %lf | ",a_epsOpt[0],a_epsOpt[1]);

  barr = 0.0;
  double fMax = -critFun(a_epsOpt, a_xx, gradObs, *a_n, *a_degree, gradf, Hess, grad, a_xGrid, *a_nv, *a_nx,barr);
  printf("f(eps)=%lf\n\n",fMax);

  ///a_epsOpt=eps;

  /// Compute Q2 @ observed data + compute empirical Mean and Variance
  double empMean[*a_degree+1];
  double empVar[*a_degree+1];
  double sumsq[*a_degree+1];
  double offSet = 0.0;


  int ii = 1;
  double q2Obs;

  for(int k=0;k<=*a_degree;k++){

      empMean[k] = 0.0;
      empVar[k] = 0.0;
      sumsq[k] = 0.0;

      double x0 = pow(a_xx[0],(double)k);

      for(int i=0;i<*a_n;i++)
      {

	  if(k==0){
	      Interp(&a_xx[i],&a_vv[i],&ii,&q2Obs,a_xGrid,a_vGrid,a_nx,a_nv,a_Q2);
	      offSet += log(q2Obs);
	  }

	  if(ABS(a_xGrid[k])>DBL_EPSILON)
	  {
	      double xi = pow(a_xx[i],(double)k);

	      empMean[k] += (xi*gradObs[i] - x0*gradObs[0]);
	      sumsq[k] += (xi*gradObs[i] - x0*gradObs[0])*(xi*gradObs[i] - x0*gradObs[0]);
	  }
      }

      empVar[k] = ((sumsq[k] - empMean[k]*empMean[k]/(*a_n))/(*a_n-1))/(*a_n);
      empMean[k] /= *a_n;
      empMean[k] += gradObs[0];
      a_normGrad[k]=empMean[k]/sqrt(empVar[k]);

    }
  offSet/=*a_n;
  *a_logL=offSet + fMax;
  
  /// Update Q2
  for(int i=0;i<*a_nv;i++)
    {
      for(int j=0;j<*a_nx;j++)
	{

	    if((a_vGrid[i] < a_xGrid[j] + *a_delta)) continue;

	    double coef=0.0;
	    if(ABS(a_xGrid[j])>DBL_EPSILON)
	    {
		for(int k=0; k<=*a_degree;k++)
		{
		    ///printf("%d ",k);
		    coef += a_epsOpt[k]*pow(a_xGrid[j],(double)k);
		}
	    }

	    a_Q2[i+*a_nv*j]*=(1.0+grad[i+*a_nv*j]*coef);
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



double critFun(double *a_eps, double *a_xx, double *a_gradObs, int a_n, int a_degree,double *a_gradf,double *a_H, double *a_gradGrid, double *a_xGrid, int a_nv, int a_nx, double a_barr){

    int deg=0;

    for(int k = 0; k<=a_degree;k++)
    {
	a_gradf[k]=0.0;
	for(int h = 0; h<=a_degree;h++)
	{
 	    a_H[k+(a_degree+1)*h]=0.0;
	}
    }
    //printf("eps: %lf %lf \n\n",a_eps[0],a_eps[1]);
    double val = 0.0;

    if(a_barr>0.0)
    {
	for(int i=0;i<a_nv;i++)
	{
	    for(int j=0;j<a_nx;j++)
	    {
		double coef=0.0;
		if(ABS(a_xGrid[j])>DBL_EPSILON)
		{
		    for(int k=0; k<=a_degree;k++)
		    {
			///printf("%d ",k);
			coef += a_eps[k]*pow(a_xGrid[j],(double)k);
		    }
		}
	    
		double denom = 1.0+a_gradGrid[i+a_nv*j]*coef;

		//if(denom<0.0) return(DBL_INF);
		if(denom<0.0) deg = 1;

		if(ABS(a_xGrid[j])>DBL_EPSILON)
		{
		    for(int k = 0; k<=a_degree;k++)
		    {
			double Dk = pow(a_xGrid[j],(double)k)*a_gradGrid[i+a_nv*j]/denom;
			a_gradf[k] -= Dk/a_barr;
		
			for(int h = 0; h<=a_degree;h++)
			{
			    a_H[k+(a_degree+1)*h] += Dk * pow(a_xGrid[j],(double)h)*a_gradGrid[i+a_nv*j]/denom/a_barr;
			}
		    }
		}
		val -= log(denom)/a_barr;

	    }
	}
    } //else a_barr=1.0;

    for(int i=0;i<a_n;i++)
    {
	double coef=0.0;
	
	for(int k = 0; k<=a_degree;k++)
	{
	    if(ABS(a_xx[i])>DBL_EPSILON)
	    {
		coef += a_eps[k]*pow(a_xx[i],(double)k);
	    }
	}

	double denom = 1.0+coef*a_gradObs[i];

	
	if(ABS(a_xx[i])>DBL_EPSILON)
	{
	    for(int k = 0; k<=a_degree;k++)
	    {
		double Dk = pow(a_xx[i],(double)k)*a_gradObs[i]/denom;
		a_gradf[k] -= Dk;///a_n;
		
		for(int h = 0; h<=a_degree;h++)
		{
		    a_H[k+(a_degree+1)*h] += Dk * pow(a_xx[i],(double)h)*a_gradObs[i]/denom;///a_n;
		}
	    }
	}
	val -= log(denom);
	//for(int m = 0; m<=a_degree;m++) val -= log(denom-a_gradMax[m]);
    }
    //val/=a_n;

    
    ///printf("Grad: %lf %lf \n\n",a_gradf[0],a_gradf[1]);

    //if(deg==1) return(DBL_INF);
    return(val);
  
}


