#include<complex>
#include <iostream>
#include <fstream>
#include <time.h> 
#include"FP.h"
#include"Global.h"
#include "mpi.h"
#include<stdlib.h>
#include<stdio.h>
#include <iomanip>
#include <gsl/gsl_linalg.h>




using namespace std;
//using Eigen::MatrixXd;
//using namespace Eigen;

FP::FP()
{
    InitialVec();
}

FP::~FP()
{
}




void FP::InitialVec()
{

//fluxS: at the half grids (i,j+1/2), and the bounday  (i,-1/2) and (i,NJs+1/2) are set to 0 
    barBS.resize(NJx+1);
    barCS.resize(NJx+1);
    barDS.resize(NJx+1);
    diffuSSDSHalfSGrid.resize(NJx+1);
    diffuSXDXHalfSGrid.resize(NJx+1);
    diffuSXHalfSGrid.resize(NJx+1);
    diffuSSHalfSGrid.resize(NJx+1);
   
    driftSHalfSGrid.resize(NJx+1);
    
    for(int i=0;i<barBS.size();i++)
    {
        barBS[i].resize(NJs+2);
        barCS[i].resize(NJs+2);
        barDS[i].resize(NJs+2);
        
        diffuSSDSHalfSGrid[i].resize(NJs+2);
        diffuSXDXHalfSGrid[i].resize(NJs+2);
        diffuSXHalfSGrid[i].resize(NJs+2);
        diffuSSHalfSGrid[i].resize(NJs+2);

        driftSHalfSGrid[i].resize(NJs+2);
    }
//...........................................................

//fluxX: at the half grids (i+1/2,j), and the bounday (-1/2,j) and (NJx+1/2,j) are set to 0 

    barBX.resize(NJx+1);
    barCX.resize(NJx+1);
    barDX.resize(NJx+1);
    
    diffuXXDX.resize(NJx+1);
    diffuXSDS.resize(NJx+1);

    
    for(int i=0;i<barBX.size();i++)
    {
        barBX[i].resize(NJs+1);
        barCX[i].resize(NJs+1);
        barDX[i].resize(NJs+1);
        
        diffuXXDX[i].resize(NJs+1);
        diffuXSDS[i].resize(NJs+1);

    }
//...................................................................


// at the  grids (i,j) drift and diffusion vector
    driftX.resize(NJx+1);
    driftS.resize(NJx+1);
    diffuXX.resize(NJx+1);
    diffuSS.resize(NJx+1);
    diffuXS.resize(NJx+1);
    diffuSX.resize(NJx+1);
    
    deltaPsiS.resize(NJx+1);
    deltaPsiX.resize(NJx+1);
    
    for(int i=0;i<driftX.size();i++)
    {
        driftX[i].resize(NJs+1);
        driftS[i].resize(NJs+1);
        
        diffuXX[i].resize(NJs+1);
        diffuSS[i].resize(NJs+1);
        diffuXS[i].resize(NJs+1);
        diffuSX[i].resize(NJs+1);
        
        deltaPsiS[i].resize(NJs+1);
        deltaPsiX[i].resize(NJs+1);
    }

}




void FP::SetFPDrift(BeamProfile &beamProfile)
{

    double deltaPOverP;

    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs;j++)
        {
            deltaPOverP = (beamProfile.js[j]-jsNominal)/jsNominal; //in the range of  [-0.2%:0.2%]*0.1

            driftX[i][j]=-0.5E0* beamProfile.jx[i];                //in the range of [0:    3E-5]*0.1  
            driftS[i][j]=-0.5E0*(beamProfile.js[j]-jsNominal);     //in the range of [-0.2%:0.2%]*0.1  
                                                                   // compare the range with Fig 3.9 and 3.10 in Fritz thesis 
        }
    }
}

void FP::SetFPDiffu(BeamProfile &beamProfile)
{

    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs;j++)
        {
            diffuXX[i][j]=1.E-10;
            diffuXS[i][j]=0.E0;
            diffuSX[i][j]=0.E0;
            diffuSS[i][j]=1.E-4;
        }
    }


}



//fluxX: at  grids (i,j), assuming two external slice (-1,j) and (Nx+1,j) are 0 to meet the bounday condition.
//the vector dim is (NJx+1,NJs+1);
void FP::SetFPBarX(BeamProfile &beamProfile,double dT)
{
 
    double dJx = (beamProfile.jxMax-beamProfile.jxMin)/NJx;

    double jsMax = (beamProfile.deltaPOverPRangeMax+1)*jsNominal;
    double jsMin = (beamProfile.deltaPOverPRangeMin+1)*jsNominal;
    double dJs=(jsMax-jsMin)/NJs;


//    dim is (NJx+1,NJs+1) 
    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs;j++)
        {
            barBX[i][j] = 0.E0;
            barCX[i][j] = 0.E0;
            barDX[i][j] = 0.E0;
            
            diffuXXDX[i][j] =0.E0;
        }
    }


    for(int i=1;i<=NJx-1;i++)
    {
        for(int j=1;j<=NJs-1;j++)
        {
            if(i==0)
            {
                diffuXXDX[i][j] = diffuXX[i+1][j] / 2/dJx;
            }
            else if(i==NJx)
            {
                diffuXXDX[i][j] = - diffuXX[i-1][j] / 2/dJx;
            }
            else
            {
                diffuXXDX[i][j] = (diffuXX[i+1][j] - diffuXX[i-1][j]) / 2/dJx;
            }
            
            
            if(j==0)
            {
                diffuXSDS[i][j] =   diffuXS[i][j+1] / 2/dJs;
            }
            else if(j==NJs)
            {
                diffuXSDS[i][j] = - diffuXS[i][j-1] / 2/dJs;
            }
            else
            {
                diffuXSDS[i][j] = (diffuXS[i][j+1] - diffuXS[i][j-1]) / 2/dJs;
            }
            

            barBX[i][j] =  driftX[i][j] - diffuXXDX[i][j]/2.E0 - diffuXSDS[i][j]/2.E0;
            barCX[i][j] = -diffuXS[i][j]/2.E0;
            barDX[i][j] = -diffuXX[i][j]/2.E0;
        }
    }
    
//    ofstream fout("testx.dat"); 
//    for(int i=0;i<=NJx;i++)
//    {
//        for(int j=0;j<=NJs;j++)
//        {
//            fout<<i<<"     "
//                <<j-NJs/2    <<"     "
//                <<barBX[i][j]<<"     "
//                <<barCX[i][j]<<"     "
//                <<barDX[i][j]<<"     "
//                <<driftX[i][j]<<"      "
//                <<diffuXSDS[i][j]<<"       "
//                <<diffuXXDX[i][j]<<"       "   // since diffuXX at the fictitious bounday are zero, it will give some values 
//                <<diffuXS[i][j]<<"     "
//                <<diffuXX[i][j]<<"     "
//                <<endl;
//        }
//    }
//    fout.close();
//    cout<<"X print"<<endl;

    int NA = (NJs+1)*(NJx+1);
    double A[NA][NA];
    
    for(int i=0;i<NA;i++)
    {
        for(int j=0;j<NA;j++)
        {
            A[i][j] = 0.E0;
        }
    }



    double coefftemp=dT/dJx;



    for(int i=2;i<=NJx-2;i++)      // (0) inner points (2:NJx-2,1:NJs-1)
    {
        for(int j=1;j<=NJs-1;j++)
        {
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] = -  barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =    barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j+1}
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =    barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] = -  barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j+1}
            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }

        }
    }
    

    
    for(int i=1;i<=1;i++)       // (1) inner ridge (1,1:NJs-1)  no i-2 terms
    {
        for(int j=1;j<=NJs-1;j++)
        {
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] = -  barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =    barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j+1}
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =    barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] = -  barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j+1}
            
        }
    }
    

    for(int i=0;i<=0;i++)       // (2) inner ridge (0,1:NJs-1) no i-2 and i-1 terms
    {
        for(int j=1;j<=NJs-1;j++)
        {
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] = -  barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =    barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j+1}
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = -  barDX[i+1][j  ]                   /(2* dJx)*coefftemp  +1;//  C_{i  ,j  }

        }
    }



    for(int i=NJx-1;i<=NJx-1;i++) //(3) inner ridge (NJX-1,1:NJs-1)  no i+2 terms
    {
        for(int j=1;j<=NJs-1;j++)
        {

            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] = -  barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =    barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j+1}
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =    barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] = -  barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j+1}
            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }

        }
    }



    for(int i=NJx;i<=NJx;i++)   //(4) inner ridge (NJX,1:NJs-1) no  i+2 and i+1 terms
    {
        for(int j=1;j<=NJs-1;j++)
        {
 
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = -  barDX[i-1][j  ]                    /(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =    barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] = -  barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j+1}
            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }
              
        }
    }
    
    for(int i=2;i<=NJx-2;i++)      //(5) inner ridge (2:NJx-2,0)  no j-1 terms
    {
        for(int j=0;j<=0;j++)
        {
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =    barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j+1}
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            

            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] = -  barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j+1}
            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }
        }
    }

    for(int i=2;i<=NJx-2;i++)      // (6) inner ridge (2:NJx-2,NJs) no j+1 terms
    {
        for(int j=NJs;j<=NJs;j++)
        {
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] = -  barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }

            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =    barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }

            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }
        }
    }
    

    
    for(int j=0;j<=0;j++)      
    {
        for(int i=0;i<=0;i++)  //  inner points (0,0) no i-2, i-1, j-1 terms 
        {
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =    barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j+1}
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ]          )        /(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
        }

        for(int i=1;i<=1;i++) //  inner points (1,0) no i-2, j-1, terms 
        {
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            

            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =    barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j+1}
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            

            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] = -  barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j+1}
        }


        for(int i=NJx-1;i<=NJx-1;i++) // inner points (NJx-1,0) no i+2, j-1, terms 
        {

            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =    barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j+1}
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] = -  barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j+1}
            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }

        }


        for(int i=NJx;i<=NJx;i++) // inner points (NJx,0) no i+2, i+1,j-1 terms 
        {

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - ( barDX[i-1][j  ]                  )/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] = -  barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j+1}
            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }

        }
    }



   for(int j=NJs;j<=NJs;j++)      
    {
        for(int i=0;i<=0;i++)  //  inner points (0,NJs) no i-2, i-1, j+1 terms 
        {
            
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] = -  barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ]                  )/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
        }

        for(int i=1;i<=1;i++) //  inner points (1,NJs) no i-2,  j+1 terms 
        {
            A[i*(NJs+1)+j][(i+2)*(NJs+1)+j  ] =    barDX[i+1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i+2,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] = -  barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }

            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =    barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }

        }

        for(int i=NJx-1;i<=NJx-1;i++) // inner points (NJx-1,NJs) no i+2, j+1 terms 
        {

            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] = -  barCX[i+1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] =    barBX[i+1][j  ]                            *coefftemp    ;//  C_{i+1,j  }
            
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - (barDX[i+1][j  ] + barDX[i-1][j  ])/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =    barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }

        }

        for(int i=NJx;i<=NJx;i++) // inner points (NJx,NJs) no i+2, i+1, j+1 terms 
        {

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = - ( barDX[i-1][j  ]                  )/(2* dJx)*coefftemp  +1;//  C_{i  ,j  }
            
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =    barCX[i-1][j  ]                   /(2* dJs)*coefftemp    ;//  C_{i-1,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] = -  barBX[i-1][j  ]                            *coefftemp    ;//  C_{i-1,j  }
            
            A[i*(NJs+1)+j][(i-2)*(NJs+1)+j  ] =    barDX[i-1][j  ]                   /(2* dJx)*coefftemp    ;//  C_{i-2,j  }
        }
    }





    double a_data[NA*NA];
    double b_data[NA];
    
    for(int i=0;i<NA;i++)
    {   
        for(int j=0;j<NA;j++)
        {
            a_data[i*NA+j] = A[i][j];
        }
    }
    
    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs;j++)
        {
            b_data[i*(NJs+1)+j] = beamProfile.psiProfile[i][j];
        }
    }
    

    gsl_matrix_view mm   = gsl_matrix_view_array (a_data, NA, NA);
    gsl_vector_view bb   = gsl_vector_view_array (b_data, NA);

    gsl_vector *x = gsl_vector_alloc (NA);

    int s;

    gsl_permutation * p = gsl_permutation_alloc (NA);
    gsl_linalg_LU_decomp (&mm.matrix, p, &s);
    gsl_linalg_LU_solve  (&mm.matrix, p, &bb.vector, x);

    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs;j++)
        {
            deltaPsiX[i][j]  =  gsl_vector_get(x,i*(NJs+1)+j) - b_data[i*(NJs+1)+j];
            beamProfile.psiProfile[i][j] = gsl_vector_get (x, i*(NJs+1)+j);
        }
    }




    gsl_permutation_free (p);
    gsl_vector_free (x);  



//    ofstream Afout("Atestx.dat"); 

//    for(int i=0;i<NA;i++)
//    {
//       for(int j=0;j<NA;j++)
//        {
//            Afout<<A[i][j]<<"     ";
//        } 
//        Afout<<endl;
//    }
//    Afout.close();
//    
//    cout<<"A_X    finsished print"<<endl;

//    getchar();
}









//fluxS: at the half grids (i,j+1/2), and the bounday  (i,-1/2) and (i,NJs+1/2) are set to 0
//the vector dim is (NJx+1,NJs+2);
void FP::SetFPBarS(BeamProfile &beamProfile,double dT)   
{
    
    double dJx = (beamProfile.jxMax-beamProfile.jxMin)/NJx;

    double jsMax = ( beamProfile.deltaPOverPRangeMax+1)*jsNominal;
    double jsMin = ( beamProfile.deltaPOverPRangeMin+1)*jsNominal;
    double dJs   = (jsMax-jsMin)/NJs;
    
    //    dim is (NJx+1,NJs+2) 
    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs+1;j++)
        {
            barBS[i][j] = 0.E0;
            barCS[i][j] = 0.E0;
            barDS[i][j] = 0.E0;
            
            driftSHalfSGrid[i][j]    =0.E0;
            
            diffuSSHalfSGrid[i][j]   =0.E0;
            diffuSSDSHalfSGrid[i][j] =0.E0;
            diffuSXDXHalfSGrid[i][j] =0.E0;
            diffuSXHalfSGrid[i][j]   =0.E0;
        }
    }



    
    for(int i=0;i<=NJx;i++)
    {
        for(int j=1;j<=NJs;j++)
        {
            if(i==0)
            {
                diffuSXDXHalfSGrid[i][j] =  (diffuSX[i+1][j-1] + diffuSX[i+1][j])/4.E0/dJx;
            }
            else if(i==NJx)
            {
                diffuSXDXHalfSGrid[i][j] = -(diffuSX[i-1][j-1] + diffuSX[i-1][j])/4.E0/dJx;
            }
            else
            {
                diffuSXDXHalfSGrid[i][j] = (diffuSX[i+1][j-1] + diffuSX[i+1][j] - diffuSX[i-1][j-1] - diffuSX[i-1][j] ) /4.E0/dJx;
            }
            

            diffuSXHalfSGrid[i][j]   = (diffuSX[i][j] + diffuSX[i][j-1]) /2.E0;
            diffuSSHalfSGrid[i][j]   = (diffuSS[i][j] + diffuSS[i][j-1])/2.E0;
            diffuSSDSHalfSGrid[i][j] = (diffuSS[i][j] - diffuSS[i][j-1])/dJs;

            driftSHalfSGrid[i][j]    = (driftS[i][j] + driftS[i][j-1])/2.E0;

            barBS[i][j] =  driftSHalfSGrid[i][j] - diffuSXDXHalfSGrid[i][j]/2.E0 - diffuSSDSHalfSGrid[i][j]/2.E0;
            barCS[i][j] = -diffuSXHalfSGrid[i][j]/2.E0;
            barDS[i][j] = -diffuSSHalfSGrid[i][j]/2.E0;
        }
        
    }

//    ofstream fout("tests.dat"); 
//    for(int i=0;i<=NJx;i++)
//    {
//        for(int j=0;j<=NJs+1;j++)
//        {
//            fout<<i<<"     "
//                <<j-NJs/2    <<"     "
//                <<barBS[i][j]<<"     "
//                <<barCS[i][j]<<"     "
//                <<barDS[i][j]<<"     "
//                <<driftSHalfSGrid[i][j]<<"      "
//                <<diffuSXHalfSGrid[i][j]<<"       "
//                <<diffuSSHalfSGrid[i][j]<<"       "
//                <<diffuSXDXHalfSGrid[i][j]<<"     "
//                <<diffuSSDSHalfSGrid[i][j]<<"     "
//                <<endl;
//        }
//    }

//    fout.close();
//    cout<<"S print"<<endl;





    int NA = (NJs+1)*(NJx+1);
    double A[NA][NA];

    for(int i=0;i<NA;i++)
    {
        for(int j=0;j<NA;j++)
        {
            A[i][j] = 0.E0;
        }
    }


    double coefftemp=dT/dJs;
//************************** solver benchmark ***************************************************************************  
//***************method-0 the longitudianl cooling with no coupling  2018-12-04--

//    for(int i=0;i<=NJx;i++)
//    {
//        for(int j=1;j<=NJs-1;j++)
//        {
//            A[i*(NJs+1)+j][(i  )*(NJs+1)+j+1] = (barBS[i  ][j+1]/2.E0 + barDS[i  ][j+1]/dJs)*coefftemp;       //  C_{i  ,j+1}
//            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
//                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1;  //  C_{i  ,j  }
//            A[i*(NJs+1)+j][(i  )*(NJs+1)+j-1] =-(barBS[i  ][j   ]/2.E0 - barDS[i  ][j  ]/dJs) *coefftemp;     //  C_{i  ,j-1}
//        }
//    }


//    for(int i=0;i<=NJx;i++)
//    {
//        for(int j=0;j<=0;j++)
//        {
//            A[i*(NJs+1)+j][(i  )*(NJs+1)+j+1] = (barBS[i  ][j+1]/2.E0 + barDS[i  ][j+1]/dJs)*coefftemp;      //  C_{i  ,j+1}
//            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
//                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }
//        }
//    }

//    for(int i=0;i<=NJx;i++)
//    {
//        for(int j=NJs;j<=NJs;j++)
//        {
//            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
//                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }
//            A[i*(NJs+1)+j][(i  )*(NJs+1)+j-1] =-(barBS[i  ][j   ]/2.E0 - barDS[i  ][j  ]/dJs) *coefftemp;    //  C_{i  ,j-1}
//        }
//    }

//***************************end of benchmark 2018-12-04 *****************************************************



//(0) get the coefficients of the ((Nx+1)*(Ns+1),(Nx+1)*(Ns+1)) matrix A. For each \Phi_{i,j}, 9 points nearby related. 

    for(int i=1;i<=NJx-1;i++)
    {
        for(int j=1;j<=NJs-1;j++)
        {
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =  barCS[i  ][j+1]/(4.E0*dJx) *coefftemp;                      //  C_{i+1,j+1}
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j+1] = (barBS[i  ][j+1]/2.E0 + barDS[i  ][j+1]/dJs)*coefftemp;      //  C_{i  ,j+1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] =- barCS[i  ][j+1]/(4.E0*dJx) *coefftemp;                      //  C_{i-1,j+1}
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] = (barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i+1,j  }
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] =-(barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i-1,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] =- barCS[i  ][j   ]/(4.E0*dJx) * coefftemp;                    //  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j-1] =-(barBS[i  ][j   ]/2.E0 - barDS[i  ][j  ]/dJs) *coefftemp;    //  C_{i  ,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =  barCS[i  ][j   ]/(4.E0*dJx) * coefftemp;                    //  C_{i-1,j-1}
        }
    }



//(1)get the coefficients of the ((Nx+1)*(Ns+1),(Nx+1)*(Ns+1)) matrix A. 

// For j=0, only Fs(i,1/2) exist;  for j=Ns, only -Fs(i,Ns-1/2) exists

    //(2.1) at the ridge (1:NJx-1,0 )   (no j-1 terms) and only Fs(i,1/2) exist;

    for(int i=1;i<=NJx-1;i++)
    {
        for(int j=0;j<=0;j++)
        {
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =  barCS[i  ][j+1]/(4.E0*dJx) *coefftemp;                      //  C_{i+1,j+1}
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j+1] = (barBS[i  ][j+1]/2.E0 + barDS[i  ][j+1]/dJs)*coefftemp;      //  C_{i  ,j+1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] =- barCS[i  ][j+1]/(4.E0*dJx) *coefftemp;                      //  C_{i-1,j+1}
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] = (barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i+1,j  }
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] =-(barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i-1,j  }
            
        }
    }




    //(2.2) at the ridge (1:NJx-1,NJS)   (no j+1 terms) and only -Fs(i,Ns-1/2) terms.   

    for(int i=1;i<=NJx-1;i++)
    {
        for(int j=NJs;j<=NJs;j++)
        {
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] = (barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i+1,j  }
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] =-(barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i-1,j  }
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] =- barCS[i  ][j   ]/(4.E0*dJx) * coefftemp;                    //  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j-1] =-(barBS[i  ][j   ]/2.E0 - barDS[i  ][j  ]/dJs) *coefftemp;    //  C_{i  ,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =  barCS[i  ][j   ]/(4.E0*dJx) * coefftemp;                    //  C_{i-1,j-1}
        }
    }





    //(2.3) at the ridge (0,1:NJs-1)  (no i-1 terms) 
    for(int i=0;i<=0;i++)
    {
      
        for(int j=1;j<=NJs-1;j++)
        {
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =  barCS[i  ][j+1]/(4.E0*dJx) *coefftemp;                      //  C_{i+1,j+1}
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j+1] = (barBS[i  ][j+1]/2.E0 + barDS[i  ][j+1]/dJs)*coefftemp;      //  C_{i  ,j+1}

            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] = (barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i+1,j  }
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }

            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] =- barCS[i  ][j   ]/(4.E0*dJx) * coefftemp;                    //  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j-1] =-(barBS[i  ][j   ]/2.E0 - barDS[i  ][j  ]/dJs) *coefftemp;    //  C_{i  ,j-1}

        }
    }
    
    //(2.4) at the ridge (NJx,1:NJs-1)  (no i+1 terms) 
    for(int i=NJx;i<=NJx;i++)
    {
        for(int j=1;j<=NJs-1;j++)
        {

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j+1] = (barBS[i  ][j+1]/2.E0 + barDS[i  ][j+1]/dJs)*coefftemp;      //  C_{i  ,j+1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] =- barCS[i  ][j+1]/(4.E0*dJx) *coefftemp;                      //  C_{i-1,j+1}
            

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] =-(barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i-1,j  }

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j-1] =-(barBS[i  ][j   ]/2.E0 - barDS[i  ][j  ]/dJs) *coefftemp;    //  C_{i  ,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =  barCS[i  ][j   ]/(4.E0*dJx) * coefftemp;                    //  C_{i-1,j-1}
        }
    }




//(3)get the coefficients of the ((Nx+1)*(Ns+1),(Nx+1)*(Ns+1)) matrix A.  at corner as \Phi_{0,0}, there are 4 points related. 
    
    for(int i=0;i<=0;i++)   
    {
        // at the coner (0,0)   only Fs(i,1/2) exist;
        for(int j=0;j<=0;j++)
        {
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j+1] =  barCS[i  ][j+1]/(4.E0*dJx) *coefftemp;                      //  C_{i+1,j+1}
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j+1] = (barBS[i  ][j+1]/2.E0 + barDS[i  ][j+1]/dJs)*coefftemp;      //  C_{i  ,j+1}
            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] = (barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i+1,j  }
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }

         }
        
        
        // at the coner (0,NJs) only Fs(i,Ns-1/2) exist;
        for(int j=NJs;j<=NJs;j++)
        {

            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j  ] = (barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i+1,j  }
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }

            
            A[i*(NJs+1)+j][(i+1)*(NJs+1)+j-1] =- barCS[i  ][j   ]/(4.E0*dJx) * coefftemp;                    //  C_{i+1,j-1}
            A[i*(NJs+1)+j][(i  )*(NJs+1)+j-1] =-(barBS[i  ][j   ]/2.E0 - barDS[i  ][j  ]/dJs) *coefftemp;    //  C_{i  ,j-1}

        }
    }
    
    
    
    for(int i=NJx;i<=NJx;i++)
    {
        // at the coner (NJx,0)   only Fs(i,1/2) exist;
        for(int j=0;j<=0;j++)
        {

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j+1] = (barBS[i  ][j+1]/2.E0 + barDS[i  ][j+1]/dJs)*coefftemp;      //  C_{i  ,j+1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j+1] =- barCS[i  ][j+1]/(4.E0*dJx) *coefftemp;                      //  C_{i-1,j+1}
            

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] =-(barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i-1,j  }

       }
        
        
        // at the coner (NJx,NJs) only Fs(i,Ns-1/2) exist;
        for(int j=NJs;j<=NJs;j++)
        {


            A[i*(NJs+1)+j][(i  )*(NJs+1)+j  ] = (barBS[i  ][j+1] - barBS[i  ][j  ])/ 2.E0     *coefftemp 
                                              - (barDS[i  ][j+1] + barDS[i  ][j  ])/dJs       *coefftemp +1; //  C_{i  ,j  }
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j  ] =-(barCS[i  ][j+1] - barCS[i  ][j  ])/(4.E0*dJx)*coefftemp;    //  C_{i-1,j  }

            A[i*(NJs+1)+j][(i  )*(NJs+1)+j-1] =-(barBS[i  ][j   ]/2.E0 - barDS[i  ][j  ]/dJs) *coefftemp;    //  C_{i  ,j-1}
            A[i*(NJs+1)+j][(i-1)*(NJs+1)+j-1] =  barCS[i  ][j   ]/(4.E0*dJx) * coefftemp;                    //  C_{i-1,j-1}

        }
    }



    double a_data[NA*NA] ;
    double b_data[NA] ;
    
    for(int i=0;i<NA;i++)
    {   
        for(int j=0;j<NA;j++)
        {
            a_data[i*NA+j] = A[i][j];
        }
    }
    
    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs;j++)
        {
            b_data[i*(NJs+1)+j] = beamProfile.psiProfile[i][j];
        }
    }
    

    gsl_matrix_view mm   = gsl_matrix_view_array (a_data, NA, NA);
    gsl_vector_view bb   = gsl_vector_view_array (b_data, NA);

    gsl_vector *x = gsl_vector_alloc (NA);

    int s;

    gsl_permutation * p = gsl_permutation_alloc (NA);
    gsl_linalg_LU_decomp (&mm.matrix, p, &s);
    gsl_linalg_LU_solve  (&mm.matrix, p, &bb.vector, x);

    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs;j++)
        {
            deltaPsiS[i][j]  =  gsl_vector_get(x,i*(NJs+1)+j) - b_data[i*(NJs+1)+j];
            beamProfile.psiProfile[i][j] = gsl_vector_get (x, i*(NJs+1)+j);
        }
    }

    gsl_permutation_free (p);
    gsl_vector_free (x);    
    
    
//    ofstream fAout("Atests.dat"); 

//    for(int i=0;i<NA;i++)
//    {
//       for(int j=0;j<NA;j++)
//        {
//            fAout<<A[i][j]<<"     ";
//        } 
//        fAout<<endl;
//    }
//    fAout.close();
//    
//    cout<<"A_S Print"<<endl;
//    getchar();
    

}




void FP::PrintData(BeamProfile beamProfile,int time)
{

    char fname[256];
    
    sprintf(fname, "beamProfilePsi_%d.dat", time);

    ofstream fout(fname);
//    fout<<" Js                      "           // 1
//        <<" Profile_Psi             "           // 2
//        <<" Drift_Psi               "           // 3
//        <<" D(Drift_Psi)/D(Js)      "           // 4
//        <<" Diffu_Psi               "           // 5
//        <<" D(Diffu_Psi)/D(Js)/D(Js)"           // 6
//        <<" D(Psi)/D(t)             "           // 7
//        <<endl;

    for(int i=0;i<=NJx;i++)
    {
        for(int j=0;j<=NJs;j++)
        {
            fout<<beamProfile.jx[i]                            <<"      "       // jx[i]                            1
                <<beamProfile.js[j]/jsNominal - 1              <<"      "       // js[j]                            2
                <<beamProfile.psiProfile[i][j]                 <<"      "       // psi[j]                           3
                <<driftX[i][j]                                 <<"      "       //
                <<driftS[i][j]                                 <<"      "       //
                <<driftX[i][j]  *beamProfile.psiProfile[i][j]  <<"      "       //
                <<driftS[i][j]  *beamProfile.psiProfile[i][j]  <<"      "       //
                <<diffuXX[i][j]                                <<"      "       //
                <<diffuXS[i][j]                                <<"      "       //
                <<diffuSX[i][j]                                <<"      "       //
                <<diffuSS[i][j]                                <<"      "       //
                <<diffuXX[i][j] *beamProfile.psiProfile[i][j]  <<"      "       //
                <<diffuXS[i][j] *beamProfile.psiProfile[i][j]  <<"      "       //
                <<diffuSX[i][j] *beamProfile.psiProfile[i][j]  <<"      "       //
                <<diffuSS[i][j] *beamProfile.psiProfile[i][j]  <<"      "       //
                <<deltaPsiS[i][j]                              <<"      "       //
                <<deltaPsiX[i][j]                              <<"      "       //
                <<endl;
        }
    }

    fout.close();
}



void FP::FPSolver(BeamProfile &beamProfile)
{

    double dT=0.01;  //seconds
    double abstime=0;

    ofstream fFPout("calculation.dat");
    for(int i=0;i<1000;i++) //       total time calculated dT*100;
    {
        beamProfile.RMSCalcu();
        
//        cout<<i<<"  "<<beamProfile.numberTest<<"     "<<beamProfile.rmsDeltapOverP<<"    "<<beamProfile.rmsEmitx<<"  "<<endl;
        
        for(int j=0; j<beamProfile.particleNumTest.size();j++)
        {
            fFPout<<beamProfile.particleNumTest[j]<<"   ";
        }
            fFPout<<endl;

        abstime = i *dT;

        SetFPDrift(beamProfile);  
        SetFPDiffu(beamProfile);
        SetFPBarS(beamProfile,dT);
        SetFPBarX(beamProfile,dT);
        
        if(i%20==0)
        {   
            PrintData(beamProfile,i);        //Print the beam profile out;
            cout<<i<<"  "<<beamProfile.numberTest<<"     "<<beamProfile.rmsDeltapOverP<<"    "<<beamProfile.rmsEmitx<<"  "<<endl;
        }

    } 

    fFPout.close();
}












