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

using namespace std; 

int mkM(vector<LONG>&m,const vector<LONG> &xs,const LONG p){    
    int degM=0;
    std::vector<LONG>linF(2,0);
    linF[1]=1;
    for(int i=0;i<xs.size();i++){
        linF[0]=(xs[i]==0?0:p-xs[i]);
        degM=pMULIP64(m,linF,degM,1,p);
    };
    return degM;
}

int main(){

    LONG p=9223372036854775783; // This is prevprime(2^63-1).
    const int degN=5;
    const int degD=5;
    const int degX=degN+degD+1;

    // Fixed size vectors for numerator and denominator.
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
    pair<vector<LONG>,int>u=newtonInterp(x,y,degX,p);
    int degU=u.second;
    // M=Product from i=0 to degX of (x-x_i).
    std::vector<LONG>m(degX+1,0);
    m[0]=1;
    int degM=mkM(m,x,p);

    printf("sizeN: %zu,sizeD: %zu,sizeX: %zu,sizeY: %zu,sizeU: %zu,sizeM: %zu\n",
    n.size(),d.size(),x.size(),y.size(),u.first.size(),m.size());

    RatReconFastWS W;
    W.init(degM+1);
    vector<LONG>rOut(W.cap,0);
    vector<LONG>tOut(W.cap,0);
    int degR=-1;
    int degT=-1;
    int flag=-999;
    const int NUM=1000000;
    auto t1=chrono::high_resolution_clock::now();

    for(int i=0;i<NUM;i++){
        flag = ratReconFastKernelWS(
        m,
        u.first,
        degX,
        degU,
        degN,
        degD,
        p,
        W,
        rOut.data(),
        degR,
        tOut.data(),
        degT
    );
      //  auto t2=chrono::high_resolution_clock::now();
      //  auto totalNano=chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
      //  double microS=totalNano/1000.0;
      //  cout<<fixed<<setprecision(3);
// cout<<"TIME (MICRO SECS.) -> "<<microS<<"\degN";
}

auto t2 = chrono::high_resolution_clock::now();

auto total_ns = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
cout << "Average ms: " << (double)total_ns / NUM / 1000.0 <<"\n"; 
   
   /* 
    auto temp2=chrono::high_resolution_clock::now();

    double avg_ns=static_cast<double>(total_ns)/1000.0;
    double avg_us=avg_ns/1000.0;
    double avg_ms=avg_ns/1.0e6;

    cout<<"TOTAL TIME: "<<total_ns<<"\degN";
    cout<<"APC (ns): "<<avg_ns<<"\degN";
    cout<<"APC (us): "<<avg_us<<"\degN";
    cout<<"APC (ms): "<<avg_ms<<"\degN";

    */

    return 0;

}
