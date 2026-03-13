#include"polyMath.h"
#include"integerMath.h"
#include"interpAlgo.h"
#include"rfrWrapper.h"
#include<vector>
#include<iostream>
#include<cstdint>
#include<time.h>
#include<chrono>

using namespace std; 

int buildModulus64(vector<LONG> &m, const vector<LONG> &xs, const LONG p){
    m.clear();
    m.push_back(1);              // m(x) = 1
    int degM = 0;

    for(int i = 0; i < (int)xs.size(); i++){
        LONG xi = xs[i] % p;
        if(xi < 0) xi += p;

        // lin(x) = x - xi = (-xi mod p) + x
        vector<LONG> lin(2);
        lin[0] = (xi == 0 ? 0 : p - xi);
        lin[1] = 1;

        degM = pMULIP64(m, lin, degM, 1, p);
        if(degM < 0) return -1;
    }

    return degM;
}

int main(){

    LONG p=9223372036854775783; // This is prevprime(2^63-1).
    const int N=5;
    const int D=5;
    const int K=N+D+1;

    vector<LONG> r(N+1,0);
    vector<LONG> t(D+1,0);
    for(int i=0;i<=N;i++){
        r[i]=rand64s(p);
        t[i]=rand64s(p);
    }

    cout<<"NUMERATOR: \n";
    dispVEC64(r);
    cout<<"\n";
    cout<<"DENOMINATOR: \n";
    dispVEC64(t);
    cout<<"\n";

    vector<LONG> x(K,0);
    for(int i=0;i<K;i++){
        x[i]=i;
    }

    vector<LONG> y(K,0);
    for(int i=0;i<K;i++){
        LONG nTemp=pEVAL64(r,N,x[i],p);
        LONG dTemp=pEVAL64(t,D,x[i],p);
        y[i]=mul64b(nTemp,modinv64b(dTemp,p),p);
    }

    cout<<"X VALUES: \n";
    dispVEC64(x);
    cout<<"\n";
    cout<<"Y VALUES: \n";
    dispVEC64(y);
    cout<<"\n";

    pair<vector<LONG>,int> u=newtonInterp(x,y,K,p);
    int degU=u.second;
    
    cout<<"U VECTOR: \n";
    dispVEC64(u.first);
    cout<<"\n";
    
    vector<LONG> m;
    int degM=buildModulus64(m,x,p);

    cout<<"M VECTOR: \n";
    dispVEC64(m);
    cout<<"\n";

    RatReconFastWS W;
    W.init(degM + degU + 4);

    vector<LONG> rOut(W.cap, 0);
    vector<LONG> tOut(W.cap, 0);

    int degR = -1;
    int degT = -1;
    int flag = -999;

    const int NUM = 100000;

auto t1 = chrono::high_resolution_clock::now();

for(int i=0; i<NUM; i++){
    flag = ratReconFastKernelWS(
        m,
        u.first,
        degM,
        degU,
        N,
        D,
        p,
        W,
        rOut.data(),
        degR,
        tOut.data(),
        degT
    );
}

auto t2 = chrono::high_resolution_clock::now();

auto total_ns = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
cout << "Average ns: " << (double)total_ns / NUM << "\n"; 
   
    
    /*
    
    auto temp2=chrono::high_resolution_clock::now();

    double avg_ns=static_cast<double>(total_ns)/1000.0;
    double avg_us=avg_ns/1000.0;
    double avg_ms=avg_ns/1.0e6;

    cout<<"TOTAL TIME: "<<total_ns<<"\n";
    cout<<"APC (ns): "<<avg_ns<<"\n";
    cout<<"APC (us): "<<avg_us<<"\n";
    cout<<"APC (ms): "<<avg_ms<<"\n";

    */

    return 0;

}
