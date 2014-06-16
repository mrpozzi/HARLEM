#define nEps 100.0
#define DBL_EPSILON 2.2204460492503131e-16

#define ABS(a)((a)>(0)?(a):(-a))
#define MIN(a,b)((a)<(b)?(a):(b))


#if !defined(M_PI)
#define M_PI 3.141592653589793238462643383280
#endif

#define DBL_EPSILON 2.2204460492503131e-16
#define DBL_INF 0x7F800000


#define bigN 1.0e7

#define alpha 0.25
#define beta 0.75
#define mu 1.0e5


void tmleStep(double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv,double *a_lambda, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv,double *a_T1,double *a_T2, double *a_Q2,double *a_normGrad, double *a_logL, double *a_epsOpt, int *a_degree);

double critFun(double *a_eps, double *a_xx, double *a_gradObs, int a_n, int a_degree,double *a_gradf,double *a_H, double *a_gradGrid, double *a_xGrid, int a_nv, int a_nx, double a_barr);

void Norm(double *a_normC, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun);
void Interp(double *a_x,double *a_y, int *a_n, double *a_Val, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun);

int LGauss(double* A,double* b,long n);
int UGauss(double* A,double* b,long n);
int Cholesky(double* A, long n);



/// Computes the normalizing constant of a function by linear interpolation
void Norm(double *a_normC, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun);

/// linear interpolation
void Interp(double *a_x,double *a_y, int *a_n, double *a_Val, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Fun);

int LGauss(double* A,double* b,long n);
int UGauss(double* A,double* b,long n);
int Cholesky(double* A, long n);




/// initialize T1 and T2
void initT(double *a_T1,double *a_T2,double *a_Stimes, double *a_Sjumps,double *a_dGn, double *a_delta, int *a_n,int *a_m, double *a_ww,  double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv);

/// Kaplan Meier
double kaplanM(double a_val,double *a_Stimes, double *a_Sjumps,int a_m);


/// initializes Q2
void initQ(double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2,int *a_N, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Q2,double *a_hOpt, double *a_Weights);
void initQslow(double *a_delta,double *a_ab, double *a_cd, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2,int *a_N, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv, double *a_Q2,double *a_hOpt, double *a_Weights);


/// Kernel Density Estimator
double kde(double a_delta,double a_tau, double a_v, double a_x, int a_n, double *a_xx, double *a_vv, double a_h[2]);
double kde_slow(double a_delta, double a_v, double a_x, int a_n, double *a_xx, double *a_vv, double a_h[2]);

/// Leave One Out approximation of the expectation of the KDE estimator
double kdeLoo(double a_delta,double a_tau, int a_n, double *a_xx, double *a_vv, double a_h[2]);
double kdeLoo_slow(double a_delta, int a_n, double *a_xx, double *a_vv, double a_h[2]);

/// Integrated Squared Error for the KDE estimator
void kdeISE(double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2,double *a_hOpt,double *a_ISE, double *a_Norm,int *a_N, double *a_Weights);
void kdeISEslow(double *a_delta,double *a_ab, double *a_cd, int *a_n, double *a_xx, double *a_vv, double *a_h1, double *a_h2, int *a_n1, int *a_n2,double *a_hOpt,int *a_N, double *a_Weights);





/// TMLE step for linear
void tmleStep1(double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv,double *a_lambda, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv,double *a_T1,double *a_T2, double *a_Q2,double *a_normGrad, double *a_logL, double *a_epsOpt);
double critFun1(double *a_xx, double *a_gradObs,int a_n,double a_eps1, double a_eps2, double a_D, double a_xD);





void tmleStep0(double *a_delta,double *a_tau, int *a_n, double *a_xx, double *a_vv,double *a_lambda, double *a_xGrid, double *a_vGrid, int *a_nx, int *a_nv,double *a_T1,double *a_T2, double *a_Q2,double *a_normGrad, double *a_logL, double *a_epsOpt);
double critFun0(double *a_xx, double *a_gradObs,int a_n,double a_eps,double a_D);
