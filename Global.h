#ifndef GLOBAL_H
#define GLOBAL_H

#include <cmath>
#include <complex>

using std::complex;

const double Epsilon 				= 8.854187817E-12;
const double C_Light 				= 299792458;                               // m/s
const double BaseCharge           =1.602176565e-19;             // C
const double BaseMass             =1.672621777e-27;             // kg
const double BaseMassIn           =931.494e+6;               // eV
const double PI                     = 4*atan(1.E0);

const  complex<double> li(0,1);
const  int NJx = 32;
const  int NJs = 16;   //meshes to for the FP solver (NJs+1)*(NJx+1); and NJx points are used for FFT. 
                      //since the approaches used, NJs>=5, NJx>=3 ate least.  
extern double jsNominal;
extern int numProcess;
extern int myRank;
#endif
