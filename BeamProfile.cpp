#include <cmath>
#include<iostream>
#include <fstream>
#include <vector>
#include "BeamProfile.h"
#include "MyFunc.h"
#include"Global.h"
#include"gsl_sf_hermite.h""

using namespace std;

BeamProfile::BeamProfile()
{
    InitialVec();
}

BeamProfile::~BeamProfile()
{
}



void BeamProfile::InitialVec()
{
    jx.resize(NJx+1);
    js.resize(NJs+1);
    
    psiProfile.resize(NJx+1);
    for(int i=0;i<psiProfile.size();i++)
    {
        psiProfile[i].resize(NJs+1);
    }

    rmsEmitxTest.resize(NJx+1);
    rmsDeltapOverPTest.resize(NJx+1);
    particleNumTest.resize(NJx+1);
}


void BeamProfile::SetInitial(LatticeRing &latticeRing)
{

    if(myRank==0)
    {
        printf("%s\n","Generate the initial beam profile...");
        cout<<endl;
    } 


    particleNum = latticeRing.particleNum;
    jsNominal = latticeRing.radius;
    jsc=(0+2.0/5.0*1.39e-3)*jsNominal; 
    jxc=7.83e-6;   

    double jxStep=(jxMax-jxMin)/NJx;
    for(int i=0;i<jx.size();i++)
    {
        jx[i]= jxMin + i * jxStep ;
    }

    double jsMax = (deltaPOverPRangeMax+1)*jsNominal;
    double jsMin = (deltaPOverPRangeMin+1)*jsNominal;
    double jsStep=(jsMax-jsMin)/NJs;
    
    for(int i=0;i<js.size();i++)
    {
        js[i]= jsMin+ i * jsStep;
    }
    
    double fnTot=0;
    for(int i=0;i<sizeof(fn)/sizeof(fn[0]);i++)
    {
        fnTot =  fnTot + fn[i];
    }


    hermiteNormIn.resize(sizeof(fn)/sizeof(fn[0]));

// do the calculation to get HermiteNormIn, which is the same as In shown below
//    for(int i=0;i<hermiteNormIn.size();i++)
//    {
//        hermiteNormIn[i] = HermiteNorm(i);
//    }

// the same as above, here is used to increae computing of the programe     
    hermiteNormIn[0]=1;
    hermiteNormIn[1]=6;
    hermiteNormIn[2]=44;
    hermiteNormIn[3]=552;
    hermiteNormIn[4]=8592;
    hermiteNormIn[5]=175200;
    hermiteNormIn[6]=4144320;




    psiNormalFactor = particleNum / PI / jsc/jxc/fnTot; 
    double q;
    for(int i=0;i<psiProfile.size();i++)
    {
        for (int j=0;j<psiProfile[i].size();j++)
        {
            q=sqrt( pow(jx[i]/jxc,2) + pow( (js[j]-jsNominal)/jsc,2) );
            psiProfile[i][j]=psiNormalFactor*PsiGen(q);
        }
    }


// here try to use the hermite_function with GSL lib, not successed yet. do it later
//    for(int i=0;i<5;i++){ 
//        cout<< gsl_sf_hermite_func(2, i)<<endl;
//        cout<< gsl_sf_hermite_func(1, i)<<endl;
//        cout<< gsl_sf_hermite_func(0, i)<<endl;
//    }
}

double BeamProfile::PsiGen(double q)
{
    double fq=0.0;
    for(int i=0;i<sizeof(fn)/sizeof(fn[0]);i++)
    {
        fq=fq +  fn[i] * HermiteFunction(2*i, q)*exp(-pow(q,2)/2)/hermiteNormIn[i];
        
//        fq=fq +  fn[i] * gsl_sf_hermite_func(2*i, q)*exp(-pow(q,2)/2)/hermiteNormIn[i];
        
    }
    
    return fq;
};

double BeamProfile::HermiteNorm(int n)
{
    double dq=1.0e-6;
    double ingerateRagne=10.;
    double temp;
    int ingerateStep = floor(ingerateRagne/dq);
    double q;
    
    for(int i=0; i<ingerateStep; i++)
    {
        q = dq * i;
        temp = temp + q * HermiteFunction(2*n, q) * exp(-pow(q,2)/2) * dq;
//         temp = temp + q * gsl_sf_hermite_func(2*n, q) * exp(-pow(q,2)/2) * dq;
    }

    return temp;
};


void BeamProfile::PrintPsi()
{

    ofstream fout("initial_profile.dat");
    fout<<"Jx"<<"  Js   "<<"deltaPOverP  "<<" Psi "<<endl;
    for(int i=0;i<psiProfile.size();i++)
    {
        for (int j=0;j<psiProfile[i].size();j++)
        {
            fout<<jx[i]<<"    "<<js[j]<<"     "<<(js[j]-jsNominal)/jsNominal<<"     "<< psiProfile[i][j]<<endl;
        }
    }
    fout.close();

}


void BeamProfile:: RMSCalcu()
{

    double emitxTot=0;
    double deltapOverPTot=0;
    double deltaJx = (jxMax-jxMin)/NJx;
    double deltaJs =  (deltaPOverPRangeMax-deltaPOverPRangeMin)/NJs;
    
    numberTest =0;
    
    for(int i=0;i<=NJx;i++)
    {
        rmsDeltapOverPTest[i] = 0.E0;
        rmsEmitxTest[i]       = 0.E0;
        particleNumTest[i]    = 0.E0;
    }
    
    
    for(int i=0;i<=NJx;i++)
    {    
        for(int j=0;j<=NJs;j++)
        {
            emitxTot        = emitxTot + jx[i] * psiProfile[i][j]*deltaJx*deltaJs;
            
            deltapOverPTot  = deltapOverPTot + pow((js[j]-jsNominal)/jsNominal,2) * psiProfile[i][j] *deltaJx*deltaJs ;
        
            numberTest = numberTest + psiProfile[i][j] *deltaJx*deltaJs;
            
            rmsDeltapOverPTest[i] = rmsDeltapOverPTest[i] + pow((js[j]-jsNominal)/jsNominal,2) * psiProfile[i][j] *deltaJx*deltaJs;
            rmsEmitxTest[i]       = rmsEmitxTest[i] + jx[i] * psiProfile[i][j]*deltaJx*deltaJs;
            particleNumTest[i]    = particleNumTest[i] + psiProfile[i][j] *deltaJx*deltaJs; 
        }
    }
    
    
    rmsEmitx        = 2*emitxTot / numberTest;
    rmsDeltapOverP  =  sqrt( deltapOverPTot/ numberTest);
    
    for(int i=0;i<=NJx;i++)
    {
        rmsDeltapOverPTest[i] = sqrt(rmsDeltapOverPTest[i] /particleNumTest[i]);
        rmsEmitxTest[i]       = 2 * rmsEmitxTest[i] / particleNumTest[i];
    }

//    cout<<rmsEmitx<<"   "<<rmsDeltapOverP<<endl;
}









