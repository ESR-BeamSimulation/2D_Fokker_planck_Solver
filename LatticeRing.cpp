#include"LatticeRing.h"
#include<iostream>
#include <fstream>
using namespace std;

LatticeRing::LatticeRing()
{
    if(myRank==0)
    {
        CodeInformaion();
    }    
    SetParam();
}

LatticeRing::~LatticeRing()
{

}


void LatticeRing::SetParam()
{
    freMid = (freLow+freHigh)/2.0;
    gammaTr = 1.0/ sqrt(mCPAlphaP);
    radius = circRing/2/PI;

    rGamma  = 1 + EK/BaseMassIn;
    rBeta   =  sqrt(1-1./pow(rGamma,2));

    ps      =  BaseMassIn/C_Light*rBeta*rGamma;

    double freq;
    freq = rBeta*C_Light/circRing;

    omegas  = 2*PI*freq;

    slipFactEta = 1.0/pow(rGamma,2) - mCPAlphaP ; 
    te      = 1.0 / 4.0 / freMid;

    harmonicH = floor(freHigh / freq);
    harmonicL = floor(freLow  / freq); 

}


void LatticeRing::CodeInformaion()
{
    cout<<"*******************************************************************"<<endl;
    cout<<"*       2D Fokker-Planck Simulation for Stochastic Cooling        *"<<endl;
    cout<<"*   This code is developed by Chao Li under the guidance from      *"<<endl;
    cout<<"*         Fritz Nolden and Mei Bai                                *"<<endl;
    cout<<"*         Any questions, please contact the authors               *"<<endl;
    cout<<"* supperli.imp@gmail.com;  f.nolden@gsi.de;   m.bai@gsi.de        *"<<endl;
    cout<<"*                 2018-11-20  FRA Univeristy-IAP                  *"<<endl;
    cout<<"*******************************************************************"<<endl;

}


