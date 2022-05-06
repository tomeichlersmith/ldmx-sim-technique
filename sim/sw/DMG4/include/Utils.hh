#pragma once

double parinv(double x, double a[], double f[], int n);


#include <iostream>

template <int NY>
double BilinearInterpolation(double X, double Y, double ArgX[], double ArgY[], double Func[][NY], int NX, int iprint)
{
  // searching for cell which bounds the point (MA, E0). Cell is numerated by its upper left index (i,j)
  int ii=-1, jj=-1;
  for (int i2 = 1; i2 < NX; i2++) {
    for (int j2 = 1; j2 < NY; j2++) {
      if( X <= ArgX[i2] && X > ArgX[i2-1] && Y <= ArgY[j2] && Y > ArgY[j2-1] ) {
        ii=i2;
        jj=j2;
        if(iprint) std::cout << "i=" << ii << " j=" << jj << " X=" << X << " Y=" << Y << std::endl;
      }
    }
  }
  if(ii < 0 || jj < 0) {std::cout << "Error in bilinear interpolation" << std::endl; exit(1);}

  // calculate interpolated value of K-factor
  double tangentX1=(X-ArgX[ii-1])/(ArgX[ii]-ArgX[ii-1]);
  double tangentY1=(Y-ArgY[jj-1])/(ArgY[jj]-ArgY[jj-1]);
  double tangentX2=(X-ArgX[ii])/(ArgX[ii-1]-ArgX[ii]);
  double tangentY2=(Y-ArgY[jj])/(ArgY[jj-1]-ArgY[jj]);
  double Result = Func[ii-1][jj-1]*tangentX2*tangentY2 +
                  Func[ii-1][jj]*tangentX2*tangentY1   +
                  Func[ii][jj-1]*tangentX1*tangentY2   +
                  Func[ii][jj]*tangentX1*tangentY1;
  //if(iprint) printf("ResultBilinear=%.4f \n", Result);
  if(iprint) std::cout << "ResultBilinear=" << Result << std::endl;
  return Result;
}
