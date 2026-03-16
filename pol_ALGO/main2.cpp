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

/*
int main(){

    LONG p=9223372036854775783; // This is prevprime(2^63-1).
    int degN=10;
    int degD=10;
    const int NUM=10000; // Total calls.
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
        std::vector<LONG> n = {7, 3, 11, 5, 2, 13, 17, 19, 23, 29, 31};
        std::vector<LONG> d = {3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

        // Populating vectors where 0<=x<p.
        
        // for(int i=0;i<=5;i++){
          //   n[i]=rand64s(p);
           //  d[i]=rand64s(p);
        // }
        // We need degN+degD+1=degX distinct points for interpolation. 
        std::vector<LONG>x(degX,0);
        for(int i=0;i<degX;i++){
            x[i]=i;
        }
        // We want f(x_i)=y_i for i=0..degX-1. 
        std::vector<LONG>y(degX,0);
        for(int i=0;i<degX;i++){
            y[i]=mul64b(pEVAL64(n.data(),degN,x[i],p),modinv64b(pEVAL64(d.data(),degD,x[i],p),p),p);
        }
        // U=Interpolation(X,Y,degX,p).
        auto interpTimeStart=chrono::steady_clock::now();
        int u=newtonInterp(x.data(),y.data(),degX,p);
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
            flag=ratReconFastKernelWS(m,y,degX,degU,
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
*/


int main(){

    LONG p=9223372036854775783; // This is prevprime(2^63-1).
    int degN=20;
    int degD=20;
    const int NUM=10000; // Total calls.
    const int iter=1;

    ofstream logFile("cppTimings.txt");
    logFile<<"Benchmark:\n";
    logFile<<"Prime p: "<<p<<"\n";
    logFile<<"Calls: "<<NUM<<"\n";
    logFile<<"Iteration degN degD InterpAvgTime avgTimeCall\n";

    for(int it=0; it<iter; it++){
        int degX=degN+degD+1;

        // Fixed size vectors for numerator and denominator.
        std::vector<LONG> n = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21
};

std::vector<LONG> d = {
    5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
};

        // Distinct interpolation points x[i] = i.
        std::vector<LONG> x(degX,0);
        for(int j=0; j<degX; j++){
            x[j]=j;
        }

        // Build original sampled values y0 once.
        std::vector<LONG> y0(degX,0);
        for(int j=0; j<degX; j++){
            LONG den = pEVAL64(d.data(),degD,x[j],p);
            if(den==0){
                std::cerr<<"Error: denominator is zero at x["<<j<<"] = "<<x[j]<<"\n";
                return 1;
            }
            LONG num = pEVAL64(n.data(),degN,x[j],p);
            y0[j]=mul64b(num,modinv64b(den,p),p);
        }

        // Warmup interpolation once.
        std::vector<LONG> y(degX,0);
        std::copy(y0.begin(), y0.end(), y.begin());
        int degU=newtonInterp(x.data(),y.data(),degX,p);

        // Time interpolation fairly over NUM calls.
        auto interpTimeStart=chrono::steady_clock::now();
        for(int k=0; k<NUM; k++){
            std::copy(y0.begin(), y0.end(), y.begin());
            degU=newtonInterp(x.data(),y.data(),degX,p);
        }
        auto interpTimeStop=chrono::steady_clock::now();
        double interpTotal=chrono::duration<double,std::micro>
            (interpTimeStop-interpTimeStart).count();
        double interpAvgTime=interpTotal/NUM;

        // Keep one interpolated copy for rat recon.
        std::vector<LONG> yInterp(degX,0);
        std::copy(y0.begin(), y0.end(), yInterp.begin());
        degU=newtonInterp(x.data(),yInterp.data(),degX,p);

        // M=Product from i=0 to degX-1 of (x-x_i).
        std::vector<LONG> m(degX+1,0);
        m[0]=1;
        int degM=mkM(m,x,p);

        RatReconFastWS W(degM);
        std::vector<LONG> rOut(degM+1,0);
        std::vector<LONG> tOut(degM+1,0);
        int degR=-1;
        int degT=-1;
        int flag=-999;

        // Warmup rat recon once.
        flag=ratReconFastKernelWS(m,yInterp,degX,degU,
            degN,degD,p,W,rOut.data(),degR,tOut.data(),degT);

        // Time rat recon over NUM calls.
        auto start=chrono::steady_clock::now();
        for(int k=0; k<NUM; k++){
            flag=ratReconFastKernelWS(m,yInterp,degX,degU,
                degN,degD,p,W,rOut.data(),degR,tOut.data(),degT);
        }
        auto stop=chrono::steady_clock::now();

        double total=chrono::duration<double,std::micro>(stop-start).count();
        double avgCallTime=total/NUM;

        logFile<<it<<"       "<<degN<<"      "
               <<degD<<"       "<<interpAvgTime<<"      "
               <<avgCallTime<<"\n";

        degN *=2;
        degD *=2;
    }

    logFile.close();
    cout<<"Complete.\n";
    return 0;
}