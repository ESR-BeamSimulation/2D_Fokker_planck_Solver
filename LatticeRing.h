#ifndef LATTICERING_H
#define LATTICERING_H

#include"Global.h"
#include <cmath>

class LatticeRing
{

public:
    LatticeRing();
    ~LatticeRing();
    
    double circRing = 108;    //meter

// fre is defined as regular frequence f 
    double freHigh = 1.7e+9;     
    double freLow  = 0.9e+9;
    double freMid;
    double te; 
    double mCPAlphaP = 0.17;
    double gammaTr;
    double slipFactEta;
    double radius;
    int harmonicH;
    int harmonicL;
    double workingQx=2.3; 
    

    
    double ZZ=92;     //Charge state
    double AA=238;    //mass number
    double EK=400.0e+6;  //Kin energy /u
    double rBeta;
    double rGamma;
    double ps;
    double omegas;
    
    
    int particleNum = 1.0e+8;

    void CodeInformaion();
    
private:
void SetParam();

};


#endif 
