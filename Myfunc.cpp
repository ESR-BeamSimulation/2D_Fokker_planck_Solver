#include "MyFunc.h"
#include <cmath>
#include <fstream>
#include<iostream>
using namespace std;



double Bessi0(double x)
{
  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i)/(tgamma(i+1))*nx;
  }
  return result;
}

double Bessi1(double x)
{
  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i+1)/(tgamma(i+1+1))*nx;
  }
  return result;
}

double Bessin(int n,double x)
{
  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i+n)/(tgamma(i+1+n))*nx;
  }
  return result;
}



double  BesselJ(int n, double x)
{

  int    m;
  double e,y;
  
  double a,b0,b1,b2;
  int i;
  a = 1;
  if (n <= 1) goto e100;
  // Calculate N! 
  for (i = 1; i < n+1; i++)  a *= i;
e100: a = 1 / a;
  if (n == 0) goto e200;
  // Calculate multiplying term
  for (i = 1; i < n+1; i++)  a *= (x / 2.0);
e200: b0 = 1.0;
  b2 = 1.0; m = 0;
  //Assemble series sum
e210: m++;
  b1 = -(x * x * b0) / (m * (m + n) * 4);
  b2 = b2 + b1;
  b0 = b1;
  // Test for convergence
  if (fabs(b1) > e)  goto e210;
  // form final answer
  y = a * b2;
  
  return y; 
}




// here to test the precise of the HermiteFunction which is good enough at least to the 12th order  
double HermiteFunction(int n, double x) 
{
    double  A[81];
    double  B[81][81];  
    int       k;
    double temp=0.0;
  int i,j;
  //Establish l0 and l1 coefficients
  B[0][0]=1.0 ; B[1][0]=0.0 ; B[1][1]=2.0;
  //Return if order is less than two
  if (n>1)  
  {
    for (i=2; i<n+1; i++)
    {
      B[i][0]=-2*(i-1)*B[i-2][0];
      for (j=1; j<i+1; j++)  
        //Basic recursion relation
        B[i][j]=2*B[i-1][j-1]-2*(i-1)*B[i-2][j];
    }
    for (i=0; i<n+1; i++)  
    {
        A[i]=B[n][i];
    
//        cout<<i<<"    "<<n<<"    "<<A[i]<<endl;
        temp =temp +  A[i] * pow(x, i);
    }
  }
  else if(n==0)
  {
        A[0]=1;
        temp =  A[0] * pow(x, 0);
  }
  else if(n==1)
  {
    A[0]=0;
    A[1]=2;
    temp = A[0] * pow(x, 0) + A[1] * pow(x, 1);
  }
  
  return temp;
}



















