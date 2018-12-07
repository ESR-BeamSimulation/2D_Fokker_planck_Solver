#include "stdio.h"
#include "string.h"
#include "mpi.h"
#include <fstream>
#include <stdlib.h> 
#include <iostream>
#include"Global.h"	
//#include"Kicker.h"
#include"BeamProfile.h"
//#include"LatticeRing.h"
//#include"TransferFunction.h"
#include"FP.h"
#include<stdlib.h> 
#include <vector>

using namespace std;
using std::vector;
using std::complex;

int main(int argc,char *argv[])
{

//    MPI_Init(&argc,&argv);
//    MPI_Status status;
//    
//    MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
//    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

    LatticeRing latticeRing;

    BeamProfile beamProfile;
    beamProfile.SetInitial(latticeRing);
    beamProfile.PrintPsi();


//    if(myRank==0)
//    {
//        beamProfile.PrintPsi(); // print the initial beam profile in Jx and Js space.
//    }

//    Kicker pickup(1,11.3,0,  6.6, 0   ,0.02, 0.07,0.025);      //mode, beta, betaPrime, disp, dispPrime, w, h, xCenter =d/2
//    Kicker kicker(1,11.3,0,  6.6, 0   ,0.02, 0.07,0.025); 
//    Kicker pickup(1,1.6, 0, 3.99, 0.   ,0.03, 0.05,0.02 );      //pickup in difference mode, kicker in sum mode
//    Kicker kicker(1,1.6, 0, 3.99, 0.   ,0.03, 0.05,0.02);     //pickup in difference mode, kicker in sum mode


//    if(myRank==0)
//    {
//        printf("Generate the kicker sensitivity as functio of x...\n");
//        pickup.PrintSensitivityToX();  //print sensitivity(x),sensitivityDx(x) as function of x
//        pickup.GetFBSenL(beamProfile,0); //print the S^L[Jx,Js],S^L_x[Jx,Js] as funciton of Jx Js;
//    }

//    TransferFunction transferFunction;
//    transferFunction.setPara(latticeRing);
//     
    FP fp;
    fp.FPSolver(beamProfile);
    
//    int m=1000;
//    int l = -1;
//    pickup.SetCorrelationCoeff(m, l, latticeRing, beamProfile);

//    fp.SetSchottkySpectrumPower(beamProfile,pickup,latticeRing);
//    fp.SetSchottkySpectrum(beamProfile,pickup,latticeRing);


//    fp.SetDrift(beamProfile, transferFunction,pickup, kicker,latticeRing);
//    fp.SetDriftDebug(beamProfile, transferFunction,pickup, kicker, latticeRing); 
//    fp.PrintData(beamProfile,2);
    
//    fp.SetDiffu(beamProfile,transferFunction,pickup, kicker, latticeRing);
    
    

//    fp.FPSolver(beamProfile, transferFunction,pickup, kicker, latticeRing);  // right now the explicitly difference method
//                                                                    // used do not ensure the conservation of particle number 
//    MPI_Finalize();

    return 0;
}











