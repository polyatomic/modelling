#ifndef _CALC_H_
#define _CALC_H_

#ifdef DECLSPEC
   #undef DECLSPEC
#endif

#ifdef _WIN32
   #ifdef MATH_EXPORTS
      #define DECLSPEC __declspec(dllexport)
   #else
      #define DECLSPEC __declspec(dllimport)
   #endif
#else
   #define DECLSPEC
#endif

#define IDX(n,i,j) (2*n-i-1)*i/2+j-i-1
#define INF 1.0e308

DECLSPEC double gaussian_kernel(double x1, double x2, double par);
DECLSPEC double laplacian_kernel(double x1, double x2, double par);
DECLSPEC double matern32_kernel(double x1, double x2, double par);
DECLSPEC double matern52_kernel(double x1, double x2, double par);
DECLSPEC double symmetrized_kernel3_0(double *x1, double *x2, double *par, double (*kernel)(double, double, double));
DECLSPEC double symmetrized_kernel3_1(double *x1, double *x2, double *par, double (*kernel)(double, double, double));
DECLSPEC double symmetrized_kernel4_0(double *x1, double *x2, double *par, double (*kernel)(double, double, double));
DECLSPEC double symmetrized_kernel4_1(double *x1, double *x2, double *par, double (*kernel)(double, double, double));
DECLSPEC double symmetrized_kernel4_2(double *x1, double *x2, double *par, double (*kernel)(double, double, double));
DECLSPEC void MatByVec(const double a[], const double b[], double c[], int n);
DECLSPEC void MatByVec(const double a[], const double b[], double c[], int m, int n);

#endif
