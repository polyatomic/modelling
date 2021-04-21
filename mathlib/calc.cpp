#include <cmath>

#ifdef _WIN32
   #include <float.h>
#endif

using namespace std;

#include "calc.h"

static double g_sqrt5 = sqrt(5.0);
static double g_5over3 = 5.0/3.0;
static double g_sqrt3 = sqrt(3.0);

double gaussian_kernel(double x1, double x2, double par)
{
   double tmp;
   tmp = x1 - x2;
   return exp(-par*tmp*tmp);;
}

double laplacian_kernel(double x1, double x2, double par)
{
   double tmp;
   tmp = fabs(x1 - x2);
   return exp(-par*tmp);;
}

double matern32_kernel(double x1, double x2, double par)
{
   double tmp;
   tmp = fabs(x1 - x2);
   return (1.0 + par*g_sqrt3*tmp)*exp(-par*g_sqrt3*tmp);
}

double matern52_kernel(double x1, double x2, double par)
{
   double tmp;
   tmp = fabs(x1 - x2);
   return (1.0 + par*g_sqrt5*tmp + par*par*g_5over3*tmp*tmp)*exp(-par*g_sqrt5*tmp);
}

double symmetrized_kernel3_0(double *x1, double *x2, double *par, double (*kernel)(double, double, double))
{
   double ret, p, k00, k11, k22, k01, k12, k21, k10, k20, k02;
   p = par[0];
   k00 = kernel(x1[0], x2[0], p);
   k11 = kernel(x1[1], x2[1], p);
   k22 = kernel(x1[2], x2[2], p);
   k01 = kernel(x1[0], x2[1], p);
   k12 = kernel(x1[1], x2[2], p);
   k21 = kernel(x1[2], x2[1], p);
   k10 = kernel(x1[1], x2[0], p);
   k20 = kernel(x1[2], x2[0], p);
   k02 = kernel(x1[0], x2[2], p);
   ret = k00*k11*k22+
         k00*k12*k21+
         k01*k10*k22+
         k01*k12*k20+
         k02*k10*k21+
         k02*k11*k20;
   return ret;
}

double symmetrized_kernel3_1(double *x1, double *x2, double *par, double (*kernel)(double, double, double))
{
   double ret, p;
   p = par[0];
   ret = (kernel(x1[0], x2[0], p)*kernel(x1[1], x2[1], p)+
          kernel(x1[0], x2[1], p)*kernel(x1[1], x2[0], p))*
         kernel(x1[2], x2[2], par[2]);
   return ret;
}

double symmetrized_kernel4_0(double *x1, double *x2, double *par, double (*kernel)(double, double, double))
{
   double ret, p, k00, k11, k22, k33, k44, k55, k12, k21, k34, k43, k01, k10, k45, k54, k20, k35, k02, k53,
          k13, k24, k31, k42, k14, k23, k32, k41, k03, k52, k04, k51, k25, k30, k15, k40, k05, k50;
   p = par[0];
   k00 = kernel(x1[0], x2[0], p);
   k11 = kernel(x1[1], x2[1], p);
   k22 = kernel(x1[2], x2[2], p);
   k33 = kernel(x1[3], x2[3], p);
   k44 = kernel(x1[4], x2[4], p);
   k55 = kernel(x1[5], x2[5], p);
   k12 = kernel(x1[1], x2[2], p);
   k21 = kernel(x1[2], x2[1], p);
   k34 = kernel(x1[3], x2[4], p);
   k43 = kernel(x1[4], x2[3], p);
   k01 = kernel(x1[0], x2[1], p);
   k10 = kernel(x1[1], x2[0], p);
   k45 = kernel(x1[4], x2[5], p);
   k54 = kernel(x1[5], x2[4], p);
   k20 = kernel(x1[2], x2[0], p);
   k35 = kernel(x1[3], x2[5], p);
   k02 = kernel(x1[0], x2[2], p);
   k53 = kernel(x1[5], x2[3], p);
   k13 = kernel(x1[1], x2[3], p);
   k24 = kernel(x1[2], x2[4], p);
   k31 = kernel(x1[3], x2[1], p);
   k42 = kernel(x1[4], x2[2], p);
   k14 = kernel(x1[1], x2[4], p);
   k23 = kernel(x1[2], x2[3], p);
   k32 = kernel(x1[3], x2[2], p);
   k41 = kernel(x1[4], x2[1], p);
   k03 = kernel(x1[0], x2[3], p);
   k52 = kernel(x1[5], x2[2], p);
   k04 = kernel(x1[0], x2[4], p);
   k51 = kernel(x1[5], x2[1], p);
   k25 = kernel(x1[2], x2[5], p);
   k30 = kernel(x1[3], x2[0], p);
   k15 = kernel(x1[1], x2[5], p);
   k40 = kernel(x1[4], x2[0], p);
   k05 = kernel(x1[0], x2[5], p);
   k50 = kernel(x1[5], x2[0], p);
   ret = k00*k11*k22*k33*k44*k55+
         k00*k12*k21*k34*k43*k55+
         k01*k10*k22*k33*k45*k54+
         k01*k12*k20*k35*k43*k54+
         k02*k10*k21*k34*k45*k53+
         k02*k11*k20*k35*k44*k53+
         k00*k13*k24*k31*k42*k55+
         k00*k14*k23*k32*k41*k55+
         k03*k10*k24*k31*k45*k52+
         k03*k14*k20*k35*k41*k52+
         k04*k10*k23*k32*k45*k51+
         k04*k13*k20*k35*k42*k51+
         k01*k13*k25*k30*k42*k54+
         k01*k15*k23*k32*k40*k54+
         k03*k11*k25*k30*k44*k52+
         k03*k15*k21*k34*k40*k52+
         k05*k11*k23*k32*k44*k50+
         k05*k13*k21*k34*k42*k50+
         k02*k14*k25*k30*k41*k53+
         k02*k15*k24*k31*k40*k53+
         k04*k12*k25*k30*k43*k51+
         k04*k15*k22*k33*k40*k51+
         k05*k12*k24*k31*k43*k50+
         k05*k14*k22*k33*k41*k50;
   return ret;
}

double symmetrized_kernel4_1(double *x1, double *x2, double *par, double (*kernel)(double, double, double))
{
   double ret, p1, p2, k00, k11, k22, k33, k44, k55, k12, k21, k34, k43, k01, k10, k45, k54, k20, k35, k02, k53;
   p1 = par[0];
   p2 = par[3];
   k00 = kernel(x1[0], x2[0], p1);
   k11 = kernel(x1[1], x2[1], p1);
   k22 = kernel(x1[2], x2[2], p1);
   k33 = kernel(x1[3], x2[3], p2);
   k44 = kernel(x1[4], x2[4], p2);
   k55 = kernel(x1[5], x2[5], p2);
   k12 = kernel(x1[1], x2[2], p1);
   k21 = kernel(x1[2], x2[1], p1);
   k34 = kernel(x1[3], x2[4], p2);
   k43 = kernel(x1[4], x2[3], p2);
   k01 = kernel(x1[0], x2[1], p1);
   k10 = kernel(x1[1], x2[0], p1);
   k45 = kernel(x1[4], x2[5], p2);
   k54 = kernel(x1[5], x2[4], p2);
   k20 = kernel(x1[2], x2[0], p1);
   k35 = kernel(x1[3], x2[5], p2);
   k02 = kernel(x1[0], x2[2], p1);
   k53 = kernel(x1[5], x2[3], p2);
   ret = k00*k11*k22*k33*k44*k55+
         k00*k12*k21*k34*k43*k55+
         k01*k10*k22*k33*k45*k54+
         k01*k12*k20*k35*k43*k54+
         k02*k10*k21*k34*k45*k53+
         k02*k11*k20*k35*k44*k53;
   return ret;
}

double symmetrized_kernel4_2(double *x1, double *x2, double *par, double (*kernel)(double, double, double))
{
   double ret, p1, p2, p3;
   p1 = par[0];
   p2 = par[1];
   p3 = par[5];
   ret = kernel(x1[0], x2[0], p1)*kernel(x1[5], x2[5], p3)*
         (kernel(x1[1], x2[1], p2)*kernel(x1[2], x2[2], p2)*kernel(x1[3], x2[3], p2)*kernel(x1[4], x2[4], p2)+
          kernel(x1[1], x2[2], p2)*kernel(x1[2], x2[1], p2)*kernel(x1[3], x2[4], p2)*kernel(x1[4], x2[3], p2)+
          kernel(x1[1], x2[3], p2)*kernel(x1[2], x2[4], p2)*kernel(x1[3], x2[1], p2)*kernel(x1[4], x2[2], p2)+
          kernel(x1[1], x2[4], p2)*kernel(x1[2], x2[3], p2)*kernel(x1[3], x2[2], p2)*kernel(x1[4], x2[1], p2));
   return ret;
}

void MatByVec(const double a[], const double b[], double c[], int n)
{
   int i, j;
   const double *v;
   for (i=0; i < n; i++)
   {
      v = a + i*n;
      c[i] = 0.0;
      for (j=0; j < n; j++)
         c[i] += v[j]*b[j];
   }
}

void MatByVec(const double a[], const double b[], double c[], int m, int n)
{
   int i, j;
   const double *v;
   for (i=0; i < m; i++)
   {
      v = a + i*n;
      c[i] = 0.0;
      for (j=0; j < n; j++)
         c[i] += v[j]*b[j];
   }
}
