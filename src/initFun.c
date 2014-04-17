#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "harlem.h"



void initT(double *a_T1,double *a_T2,double *a_Stimes, double *a_Sjumps,double *a_dGn, double *a_delta, int *a_n,int *a_m, double *a_ww,  double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv)
{
    for(int i=0; i<*a_nv; i++)
    {
        for(int j=0; j<*a_nx; j++)
        {
	    a_T1[i + *a_nv * j] = 0.0;
	    a_T2[i + *a_nv * j] = 0.0;

            if(a_vGrid[i] > a_xGrid[j] + *a_delta)
            {
                double denominatore = 0.0;
                int k = 0;
		
                while((a_vGrid[i] > a_ww[k]) && (k < *a_n))
                {
                    if(a_ww[k] > a_xGrid[j]) denominatore += a_dGn[k] * kaplanM(a_vGrid[i] - a_ww[k], a_Stimes, a_Sjumps, *a_m);
                    k++;
                }
                if(denominatore > sqrt(DBL_EPSILON))
                {
                    a_T1[i + *a_nv * j] = (a_vGrid[i] - a_xGrid[j])/denominatore;
                    a_T2[i + *a_nv * j] = 1/denominatore;
                }
            } 
        }
    }
}

double kaplanM(double a_val, double *a_Stimes, double *a_Sjumps, int a_m)
{
    if(a_val<a_Stimes[0]) return(1.0);
    if(a_val>a_Stimes[a_m-1]) return(a_Sjumps[a_m-1]);
    
    int up = a_m - 1;
    int low = 0;
    int i = ceil((up+low)/2);
    while(1)
    {
        
        if(a_val == a_Stimes[i+1]) return(a_Sjumps[i]);
        if(a_val == a_Stimes[i]) return(a_Sjumps[i-1]);
        if(a_val == a_Stimes[i-1]) return(a_Sjumps[i-1]);
        
        if(a_Stimes[i]>a_val)
        {
            if(a_Stimes[i-1] <= a_val)
            {
                return(a_Sjumps[i-1]);
            }else
            {
                up = i;
                i = floor((up+low)/2);
            }
        }else if(a_Stimes[i+1]>a_val)
        {
            return(a_Sjumps[i]);
        }else
        {
            low = i;
            i = ceil((up+low)/2);
        }
        
    }
}





void initQ(double *a_delta,double *a_ab, double *a_cd, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2,int *a_N, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Q2,double *a_hOpt, double *a_Weights)
{
    if(*a_n1>1 || *a_n2>1)
    {
	kdeISE(a_delta, a_ab, a_cd, a_n, a_xx, a_vv, a_h1, a_h2, a_n1, a_n2, a_hOpt, a_N, a_Weights);
    } else {
	a_hOpt[0] = *a_h1;
	a_hOpt[1] = *a_h2;
    }
    for(int i=0; i<*a_nv;i++)
    {
        for(int j=0;j<*a_nx;j++)
        {
	    if(i==0 || j==0)
	    {
		a_Q2[i + *a_nv * j] = 0.0;
		continue;
	    }
	    double y = (a_vGrid[i] - a_xGrid[j] - *a_delta);
	    if(y <= 0)
	    {
		a_Q2[i + *a_nv * j]=0.0;
		continue;
	    }
	    a_Q2[i + *a_nv * j] = kde(*a_delta, y, a_xGrid[j], *a_n, a_xx, a_vv, a_hOpt) / (a_xGrid[j] * y);
        }
    }

}


void kdeISE(double *a_delta, double *a_ab, double *a_cd, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2,double *a_hOpt, int *a_N, double *a_Weights)
{
    
    double ISE = 0.0;
    double min_ISE = 0.0;

    
    for (int i=0; i<*a_n1; i++)
    {
        
        double h[2];
        h[1]=a_h2[i];
	double c = a_cd[0] - h[1];
	double d = a_cd[1] + h[1];
	double stepU2 = (d - c)/((double)*a_N);
        
        for (int j=0; j<*a_n2; j++)
        {
            
            double avgSq = 0.0;
            h[0]=a_h1[j];

	    printf("%lf %lf\n",h[0],h[1]);

	    double a = a_ab[0] - h[0];
	    double b = a_ab[1] + h[0];
	    double stepU1 = (b - a)/((double)*a_N);
            
            /// Simpson Formula in 2D
	    for(int l=0;l<=*a_N;l++)
            {
		for(int k=0;k<=*a_N;k++)
                {

		    double term1 = kde(*a_delta,c+l*stepU2,a+k*stepU1,*a_n,a_xx,a_vv, h);
		    //double term1 = kde(*a_delta,*a_tau,c+l*step-(a+k*step)- *a_delta,a+k*step,*a_n,a_xx,a_vv, h);
		    avgSq += a_Weights[l+(*a_N+1)*k]*term1*term1;
		    //kde(double a_delta,double a_tau, double a_v, double a_x, int a_n, double *a_xx, double *a_vv, double a_h[2])
		}
            }

            avgSq *= stepU1*stepU2/9.0;
            
            double loo=kdeLoo(*a_delta,*a_n,a_xx,a_vv,h);
                        
            ISE = avgSq-2*loo;

	    
            if(i==0 && j==0)
            {
                min_ISE = ISE;
                a_hOpt[0]=h[0]; a_hOpt[1]=h[1];
            }else if (ISE < min_ISE) 
            {
                min_ISE = ISE;
                a_hOpt[0]=h[0]; a_hOpt[1]=h[1];
            }
            
        }
    }
    
}


double kde(double a_delta, double a_v, double a_x, int a_n, double *a_xx, double *a_vv, double a_h[2])
{
  
  if(a_v <= 0.0) return(0.0);
  if(a_x <= 0.0) return(0.0);
  double x = log(a_x); double y = log(a_v);
  
  double p = 0.0;
  
  for (int i=0; i<a_n; i++) 
  {
   
    double yy = a_vv[i] - a_xx[i] - a_delta;

    if(yy <= 0.0) continue;
    if(a_xx[i] <= 0.0) continue;

    double z = (x - log(a_xx[i]))/a_h[0];
    double t = (y - log(yy))/a_h[1];
    

    if((t<1.0) && (t>-1.0) && (z<1.0) && (z>-1.0)) p += (1-(t*t)) * (1-(z*z)); 

  }

  p /= a_n * (a_h[0] * a_h[1]) * 16/9;
  
  return(p);
	
}


double kdeLoo(double a_delta, int a_n, double *a_xx, double *a_vv, double a_h[2])
{
  
  double LOO=0.0;
  
  for (int i=0; i<a_n; i++) 
  {
  
    double y = a_vv[i] - a_xx[i] - a_delta;
    double x = a_xx[i];
    
    if(y <= 0.0) continue;
    if(x <= 0.0) continue;

    x = log(x); y = log(y);
    
    for (int j=0; j<a_n; j++) 
    {
      
      if(i==j)continue;

      double yy = a_vv[j]-a_xx[j]-a_delta;

      if(yy <= 0.0) continue;
      if(a_xx[j] <= 0.0) continue;

      double z = (x-log(a_xx[j]))/a_h[0];
      double t = (y-log(yy))/a_h[1];

      if((t<1.0)&&(t>-1.0)&&(z<1.0)&&(z>-1.0)) LOO+=(1-(t*t))*(1-(z*z));

    }

  }
	
  LOO/=a_n*(a_n-1)*(a_h[0]*a_h[1])*16/9;

  return(LOO);
  	
}
