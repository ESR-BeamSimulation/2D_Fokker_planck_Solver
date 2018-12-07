
#ifndef MYFUNC_H
#define MYFUNC_H

#include <vector>
#include <string>
#include <ostream>


using std::string;
using std::vector;


//void   StringSplit(const string &s,vector<string> &vec);
//int    GetLineNumber(const string &p);

//first modified Bessel function Bessin(int n, double x).
double Bessi0(double x);
double Bessi1(double x);
double Bessin(int n,double x);


//first kind of bessle function  BesselJ (int n, double x)
double BesselJ(int n, double x);

double HermiteFunction(int n, double x); 



//void PrintV1(std::ostream &par,const        vector<double>    &particle);
//void PrintV2(std::ostream &par,const vector<vector<double>  > &particle);
#endif
