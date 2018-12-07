#ifndef BEAMPROFILE_H
#define BEAMPROFILE_H

#include "Global.h"
#include <vector>
#include <complex>

using namespace std;
using std::vector;
using std::complex;

#include "LatticeRing.h"

class BeamProfile
{


public:
    BeamProfile();
    ~BeamProfile();



    double psiNormalFactor;


    double jsc;  
    double jxc;

    double numberTest;
    
    double jxMax= 3.0e-5;
    double jxMin=-0.0e-5;
    double deltaPOverPRangeMax=2E-3;
    double deltaPOverPRangeMin=-2E-3;

    double fn[7]={1.3949,1.7827,1.0047,0.2033,-0.1729,-0.0209,0.0536};
    double particleNum;
    vector<double > jx;
    vector<double > js;
    vector<vector<double> > psiProfile;
    vector<double > hermiteNormIn;
    double PsiGen(double q);
    double HermiteNorm(int n);
    
    double rmsEmitx;
    double rmsDeltapOverP;
    
    vector <double > rmsEmitxTest;
    vector <double > rmsDeltapOverPTest;
    vector <double > particleNumTest;
    
    void InitialVec();

    void RMSCalcu();
    void SetInitial(LatticeRing &latticeRing);
    void PrintPsi();

private:

};



#endif
