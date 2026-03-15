#include"polyMath.h"
#include"integerMath.h"
#include"interpAlgo.h"
#include"rfrWrapper.h"
#include<vector>
#include<iostream>
#include<cstdint>
#include<time.h>
#include<chrono>
#include<iomanip>
#include<fstream>

using namespace std; 

int mkM(vector<LONG>&m,const vector<LONG> &xs,const LONG p){    
    int degM=0;
    std::vector<LONG>linF(2,0);
    linF[1]=1;
    for(int i=0;i<xs.size();i++){
        linF[0]=(xs[i]==0?0:p-xs[i]);
        degM=pMULIP64(m.data(),linF.data(),degM,1,p);
    };
    return degM;
}

int main(){

    LONG p=9223372036854775783; // This is prevprime(2^63-1).
    int degN=160;
    int degD=160;
    const int NUM=1000000; // Total calls.
    const int iter=1;
    ofstream logFile("cppTimings.txt");
    logFile<<"Benchmark:\n";
    logFile<<"Prime p: "<<p<<"\n";
    logFile<<"Calls: "<<NUM<<"\n";
    logFile<<"Iteration degN degD InterpTime avgTimeCall\n";

    for(int i=0;i<iter;i++){
        int degX=degN+degD+1;
        // Fixed size vectors for numerator and denominator.
        // auto vecBuildTimeStart=chrono::steady_clock::now();
        std::vector<LONG>n(degN+1,0);
        std::vector<LONG>d(degD+1,0);
        // Populating vectors where 0<=x<p.
        for(int i=0;i<=degN;i++){
            n[i]=rand64s(p);
            d[i]=rand64s(p);
        }
        // We need degN+degD+1=degX distinct points for interpolation. 
        std::vector<LONG>x(degX,0);
        for(int i=0;i<degX;i++){
            x[i]=i;
        }
        // We want f(x_i)=y_i for i=0..degX-1. 
        std::vector<LONG>y(degX,0);
        for(int i=0;i<degX;i++){
            y[i]=mul64b(pEVAL64(n,degN,x[i],p),modinv64b(pEVAL64(d,degD,x[i],p),p),p);
        }
        // U=Interpolation(X,Y,degX,p).
        vector<LONG> b(degX);
        auto interpTimeStart=chrono::steady_clock::now();
        int u=newtonInterp(x.data(),y.data(),degX,p,b.data());
        auto interpTimeStop=chrono::steady_clock::now();
        double interpTimeFinal=chrono::duration<double,
        std::micro>(interpTimeStop-interpTimeStart).count();
        int degU=u;
        // M=Product from i=0 to degX of (x-x_i).
        std::vector<LONG>m(degX+1,0);
        m[0]=1;
        int degM=mkM(m,x,p);
        RatReconFastWS W(degM);
        vector<LONG>rOut(degM+1,0);
        vector<LONG>tOut(degM+1,0);
        int degR=-1;
        int degT=-1;
        int flag=-999;
        // auto vecBuildTimeStop=chrono::steady_clock::now();
        // double vecTotalTime=chrono::duration<double,
        // std::micro>(vecBuildTimeStop-vecBuildTimeStart).count();
        auto start=chrono::steady_clock::now();
        for(int i=0;i<NUM;i++){
            flag=ratReconFastKernelWS(m,b,degX,degU,
            degN,degD,p,W,rOut.data(),degR,tOut.data(),
            degT);
        }
        auto stop=chrono::steady_clock::now();
        auto total=chrono::duration_cast<chrono::microseconds>
        (stop-start).count();
        double avgCallTime=static_cast<double>(total)/NUM;

        logFile<<i<<"       "<<degN<<"      "
        <<degD<<"       "<<interpTimeFinal<<"      "
        <<avgCallTime<<"\n";
        degN *=2;
        degD *=2;
    }
    logFile.close();
    cout<<"Complete.\n";
    return 0;
}
