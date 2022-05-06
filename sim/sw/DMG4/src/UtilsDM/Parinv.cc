/* +++++++++++++++++++ Util routines: Spline interpolation ++++++++++++ */
/* ++++++++++++++++++++++ C Include Files ++++++++++++++++++++++ */
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define EPSPARINV 1.e-8

double parinv(double x, double a[], double f[], int n)
{
//
//    Interpolation at the point x. Function f(a) is tabulated
//    in arrays a, f with dimension n.
//
  int k1, k2, k3;

  if(n < 3) {std::cerr << "parinv: insufficient number of points" << std::endl; exit(1);}
  if(x < a[0]) {
    double c = fabs(x - a[0]);
    if(c < EPSPARINV*fabs(a[1]-a[0])) return a[0];
    k1 = 0;
  }
  else if(x > a[n-1]) {
    double c = fabs(x - a[n-1]);
    if(c < EPSPARINV*fabs(a[n-1]-a[n-2])) return a[n-1];
    k1 = n-3;
  }
  else {
    k1 = 0;
    k2 = n-1;
    k3 = k2 - k1;
    while(k3 > 1) {
      k3 = k1 + k3/2;
      if( a[k3]-x == 0 ) return f[k3];
      if( a[k3]-x < 0 ) k1 = k3;
      if( a[k3]-x > 0 ) k2 = k3;
      k3 = k2 - k1;
    }
    if(k2 == n-1) k1 = n - 3;
  }
  if(k1 < 0 || k1 > n-3) {std::cerr << "parinv: wrong index found" << std::endl; exit(1);}
  double b1 = a[k1];
  double b2 = a[k1+1];
  double b3 = a[k1+2];
  double b4 = f[k1];
  double b5 = f[k1+1];
  double b6 = f[k1+2];
  return b4 * ((x-b2)*(x-b3))/((b1-b2)*(b1-b3)) +
         b5 * ((x-b1)*(x-b3))/((b2-b1)*(b2-b3)) +
         b6 * ((x-b1)*(x-b2))/((b3-b1)*(b3-b2));
}
