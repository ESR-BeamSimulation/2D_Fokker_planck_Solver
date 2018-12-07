#ifndef FP_H
#define FP_H
#include"BeamProfile.h"
//#include"Kicker.h"
//#include"TransferFunction.h"
#include"LatticeRing.h"
#include <complex>
#include <vector>


using namespace std;
using std::vector;
using std::complex;

class FP
{
public:
    FP();
    ~FP();

//    void  SetDrift(BeamProfile &beamProfile, TransferFunction &transferFunction,Kicker &pickup, Kicker &kicker, LatticeRing &latticeRing);
//    void  SetDriftDebug(BeamProfile &beamProfile, TransferFunction &transferFunction,Kicker &pickup, Kicker &kicker, LatticeRing &latticeRing);
//    
//    void  SetDiffu(BeamProfile &beamProfile, TransferFunction &transferFunction,Kicker &pickup, Kicker &kicker, LatticeRing &latticeRing);
//    void  SetDiffuDebug(BeamProfile &beamProfile, TransferFunction &transferFunction,Kicker &pickup, Kicker &kicker, LatticeRing &latticeRing);
//    
//    void FPSolver(BeamProfile &beamProfile, TransferFunction &transferFunction,Kicker &pickup, Kicker &kicker, LatticeRing &latticeRing);

    vector<vector<double > > driftX;
    vector<vector<double > > driftS;
    vector<vector<double > > diffuXX;
    vector<vector<double > > diffuSS;
    vector<vector<double > > diffuXS;
    vector<vector<double > > diffuSX;
    
// used in the fluxX calculation
//    vector<vector<double > > diffuXXDXHalfXGrid;
//    vector<vector<double > > diffuXSDSHalfXGrid;
//    vector<vector<double > > diffuXXHalfXGrid;
//    vector<vector<double > > diffuXSHalfXGrid;

//    vector<vector<double > > driftXHalfXGrid;

    vector<vector<double > > barBX;
    vector<vector<double > > barCX;
    vector<vector<double > > barDX;
    
    vector<vector<double > > diffuXXDX;
    vector<vector<double > > diffuXSDS;


// used in the fluxS calculation
    vector<vector<double > > diffuSSDSHalfSGrid;
    vector<vector<double > > diffuSXDXHalfSGrid;
    vector<vector<double > > diffuSXHalfSGrid;
    vector<vector<double > > diffuSSHalfSGrid;

    vector<vector<double > > driftSHalfSGrid;

    vector<vector<double > > barBS;
    vector<vector<double > > barCS;
    vector<vector<double > > barDS;


    vector<vector<double > > deltaPsiS;
    vector<vector<double > > deltaPsiX;
    


//    void SetSchottkySpectrum(BeamProfile &beamProfile,Kicker &pickup,LatticeRing &latticeRing);
//    void SetSchottkySpectrumPower(BeamProfile &beamProfile,Kicker &pickup,LatticeRing &latticeRing);

      void InitialVec();


    void SetFPDrift(BeamProfile &beamProfile);
    void SetFPDiffu(BeamProfile &beamProfile);
    void SetFPBarX(BeamProfile &beamProfile,double dT);
    void SetFPBarS(BeamProfile &beamProfile,double dT);
    void FPSolver(BeamProfile &beamProfile);


private:
//    void TimeDisplay();

    void PrintData(BeamProfile beamProfile,int time);
    void SetFPDletaPsiS();
    void SetFPDletaPsiX();
//    void Push(BeamProfile &beamProfile, double &timeStep);

};




#endif
