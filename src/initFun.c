#include <stdlib.h> 
#include <stdio.h> 
#include<math.h>
///#include <time.h> 

#if !defined(M_PI)
#  define M_PI 3.141592653589793238462643383280
#endif

#define DBL_EPSILON 2.2204460492503131e-16


void initT(double *a_T1,double *a_T2,double *a_Stimes, double *a_Sjumps,double *a_dGn, double *a_delta, int *a_n,int *a_m, double *a_ww,  double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv);
double kaplanM(double a_val,double *a_Stimes, double *a_Sjumps,int a_m);



void initQ(int  *a_Order, double *a_Nodes, double *a_Weights,double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Q2,double *a_hOpt);
double kde(double a_delta,double a_tau, double a_x, double a_v, int a_n, double *a_xx, double *a_vv, double a_h[2]);
double kdeLoo(double a_delta,double a_tau, int a_n, double *a_xx, double *a_vv, double a_h[2]);
void kdeISE(int  *a_Order, double *a_Nodes, double *a_Weights,double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2,double *a_hOpt,double *a_ISE, double *a_Norm);


void initT(double *a_T1,double *a_T2,double *a_Stimes, double *a_Sjumps,double *a_dGn, double *a_delta, int *a_n,int *a_m, double *a_ww,  double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv)
{

  for(int i=0; i<*a_nv;i++)
    {
      for(int j=0;j<*a_nx;j++)
	{
	  //printf("%d %d\n",i,j);
	  if(a_vGrid[i]>= a_xGrid[j]+*a_delta)
	    {
	      double denominatore = 0.0;
	      int k = 0;
	      while((a_vGrid[i]>a_ww[k])&&(k<*a_n))
		{
		  if(a_ww[k]>a_xGrid[j]) denominatore += a_dGn[k]*kaplanM(a_vGrid[i]-a_ww[k],a_Stimes, a_Sjumps,*a_m);
		  k++;
		}
	      if(denominatore<=sqrt(DBL_EPSILON))
		{
		  a_T1[i+*a_nv*j]=0.0;
		  a_T2[i+*a_nv*j]=0.0;
		}else
		{
		  a_T1[i+*a_nv*j]= (a_vGrid[i]- a_xGrid[j])/denominatore;
		  a_T2[i+*a_nv*j]=1/denominatore;
		}
	    }else
	    {
	      a_T1[i+*a_nv*j]=0.0;
	      a_T2[i+*a_nv*j]=0.0;
	    }
	}
    }
}

double kaplanM(double a_val,double *a_Stimes, double *a_Sjumps,int a_m)
{
  if(a_val<a_Stimes[0]) return(1.0);
  if(a_val>a_Stimes[a_m-1]) return(a_Sjumps[a_m-1]);
  //int found=0;

  int up=a_m-1;
  int low=0;
  int i = ceil((up+low)/2);
  while(1)
    {

      //if(a_Stimes[i]<a_val&&a_Stimes[i+1]>=a_val)return(a_Sjumps[i]);
      //i++;
      
      if(a_val==a_Stimes[i+1])return(a_Sjumps[i]);
      if(a_val==a_Stimes[i])return(a_Sjumps[i-1]);
      if(a_val==a_Stimes[i-1])return(a_Sjumps[i-1]);
      
      if(a_Stimes[i]>a_val)
	  {
	    if(a_Stimes[i-1]<=a_val)
	    {
	      return(a_Sjumps[i-1]);
	    }else 
	    {
	      up=i;
	      i = floor((up+low)/2);
	    }
       }else if(a_Stimes[i+1]>a_val)
       {
	 return(a_Sjumps[i]);
       }else 
       {
	 low=i;
	 i = ceil((up+low)/2);
       }
      
    }
      
}

void initQ(int  *a_Order, double *a_Nodes, double *a_Weights,double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Q2,double *a_hOpt)
{

  double ISE=0.0;
  //double hOpt[2];
  double normC=0.0;
  if(*a_n1>1||*a_n2>1)
    {
      kdeISE(a_Order, a_Nodes, a_Weights,a_delta,a_tau,a_n,a_xx, a_vv, a_h1, a_h2, a_n1, a_n2,a_hOpt,&ISE,&normC);
    }else
    {

      double coeff=(*a_tau-*a_delta)/2.0;
      double shift=(*a_tau+*a_delta)/2.0;
	 
      a_hOpt[0]=a_h1[0];
      a_hOpt[1]=a_h2[0];


      for(int l=0;l<*a_Order;l++){
	    
	double y[2];
	y[0] = coeff*a_Nodes[l]+shift;
	
	for(int k=0;k<*a_Order;k++){
	  
	  y[1] = coeff*a_Nodes[k]+shift;
	  
	  double fHat=kde(*a_delta,*a_tau,y[0],y[1],*a_n,a_xx,a_vv,a_hOpt);
	  
	  normC+=a_Weights[l]*a_Weights[k]*fHat;
	  
	}


      }

      normC*=coeff*coeff;

    }
   
  //printf("C=%g\n",normC);

  for(int i=0; i<*a_nv;i++)
    {
      for(int j=0;j<*a_nx;j++)
	{
	  a_Q2[i+*a_nv*j]=kde(*a_delta, *a_tau,a_vGrid[i]- a_xGrid[j], a_vGrid[i], *a_n,a_xx,a_vv, a_hOpt)/normC;
	}
    }
}


void kdeISE(int  *a_Order, double *a_Nodes, double *a_Weights,double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2,double *a_hOpt,double *a_ISE, double *a_Norm)
{

  //srand((unsigned)time(NULL));
  //double coeff=(*a_tau-*a_delta);
  //double shift=*a_delta;

  double coeff=(*a_tau-*a_delta)/2.0;
  double shift=(*a_tau+*a_delta)/2.0;
  double ISE=0.0;

  for (int i=0; i<*a_n1; i++) {

    double h[2];
    h[0]=a_h1[i];

    for (int j=0; j<*a_n2; j++) {
            
      double avgSq = 0.0;
      double normC = 0.0;

      h[1]=a_h2[j];

      for(int l=0;l<*a_Order;l++){
	
	double y[2];
	y[0] = coeff*a_Nodes[l]+shift;

	for(int k=0;k<*a_Order;k++){

	  y[1] = coeff*a_Nodes[k]+shift;
	 
	  double fHat=kde(*a_delta,*a_tau,y[0],y[1],*a_n,a_xx,a_vv, h);

	  normC+=a_Weights[l]*a_Weights[k]*fHat;
	  avgSq+=a_Weights[l]*a_Weights[k]*fHat*fHat;

	}
      }

      normC*=coeff*coeff;

      avgSq*=coeff*coeff;
      avgSq/=normC*normC;
      
      double loo=kdeLoo(*a_delta,*a_tau,*a_n,a_xx,a_vv,h);
      
      loo/=normC;

      ISE = avgSq-2*loo;

      //printf("[h1=%lf; h2=%lf] -> %g - %g = ISE=%g (%g) \n",a_h1[i],a_h2[j],avgSq,2*loo,ISE,*a_ISE);
      

      if(i==0&&j==0){
	//printf("FOUND!! fsq=%g; loo=%g C=%g\n",avgSq,loo,normC);
	*a_ISE=ISE;
	a_hOpt[0]=h[0]; a_hOpt[1]=h[1];
	*a_Norm=normC;
      }else if (ISE<*a_ISE) {
	//printf("FOUND!! fsq=%g; loo=%g C=%g\n",avgSq,loo,normC);
	*a_ISE=ISE;
	a_hOpt[0]=h[0]; a_hOpt[1]=h[1];
	*a_Norm=normC;
      }
      
    }
  }
  
  //printf("%g\n",*a_Norm);	
  
}


double kde(double a_delta,double a_tau, double a_x, double a_v, int a_n, double *a_xx, double *a_vv, double a_h[2])
{
  
  //double y=a_v-a_x;
  double y=a_x;
  
  if((y<a_delta)||(y>a_tau)) return(0);
  if((a_v<a_delta)||(a_v>a_tau)) return(0);
  
  double p=0.0;
  
  for (int i=0; i<a_n; i++) {
   
    double yy=a_vv[i]-a_xx[i];
    double t=(y-yy)/a_h[0];
    double z=(a_v-a_vv[i])/a_h[1];

    //p+=exp(-(t*t)/2-(z*z)/2);
    if((t<1.0)&&(t>-1.0)&&(z<1.0)&&(z>-1.0)) p+=(1-(t*t))*(1-(z*z)); 

  }

  //p/=a_n*a_h[0]*a_h[1]*(2*M_PI);
  p/=a_n*(a_h[0]*a_h[1])*16/9;
  
  return(p);
	
}


double kdeLoo(double a_delta,double a_tau, int a_n, double *a_xx, double *a_vv, double a_h[2])
{
  
  double LOO=0.0;
  
  for (int i=0; i<a_n; i++) {
  
    //double loo=0.0;

    //double y=a_xx[i];
    double y=a_vv[i]-a_xx[i];
    double v=a_vv[i];
    
    if((y<a_delta)||(y>a_tau)) continue;
    if((v<a_delta)||(v>a_tau)) continue;
    
    for (int j=0; j<a_n; j++) {
      
      if(i==j)continue;

      //double yy=a_xx[j];
      double yy=a_vv[j]-a_xx[j];
      double t=(y-yy)/a_h[0];
      double z=(v-a_vv[j])/a_h[1];
      //LOO+=exp(-(t*t)/2-(z*z)/2);
      if((t<1.0)&&(t>-1.0)&&(z<1.0)&&(z>-1.0)) LOO+=(1-(t*t))*(1-(z*z));

    }
    //loo/=(a_n-1)*a_h[0]*a_h[1]*(2*M_PI);
    //loo/=(a_n-1)*(a_h[0]*a_h[1])*16/9;
    //LOO+=loo;
  }
	
  LOO/=a_n*(a_n-1)*(a_h[0]*a_h[1])*16/9;
  //LOO/=a_n*(a_n-1)*(a_h[0]*a_h[1])*(2*M_PI);
  return(LOO);
  	
}
