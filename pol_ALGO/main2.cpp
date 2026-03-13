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

    RatReconFastWS W(degM);
    vector<LONG>rOut(degM+1,0);
    vector<LONG>tOut(degM+1,0);
    int degR=-1;
    int degT=-1;
    int flag=-999;
    const int NUM=100000;
    auto start=chrono::high_resolution_clock::now();

    for(int i=0;i<NUM;i++){
        flag=ratReconFastKernelWS(m,u.first,degX,degU,
        degN,degD,p,W,rOut.data(),degR,tOut.data(),
        degT);
    }
    printf("%d\n",flag);
    auto stop=chrono::high_resolution_clock::now();
    auto total=chrono::duration_cast<chrono::microseconds>
    (stop-start).count();
    std::cout<<"Average micro(s): "<<(double)total/NUM<<"\n"; 
    return 0;
}
