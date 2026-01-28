#include<bits/stdc++.h> 
#include"integerMath.h"
#include"helperF.h"
#include"polyMath.h"

using namespace std; 

using LONG=int_fast64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;
using clockk=chrono::steady_clock;

static volatile LONG sink=0;

template<class F>
double benchMark(const char*name,int iters,F&& fn){
//    for(int i=0;i<max(1,iters/10);++i) fn();
    auto t0=clockk::now();
    for(int i=0;i<iters;++i) fn();
    auto t1=clockk::now();
    double us=chrono::duration<double,micro>(t1-t0).count();
    cout<<name<< " -> " << (us/iters) << " us per ops.\n";
    return us/iters;
}

int main(){
    const LONG p=((1ULL<<61)-1); // This is 2^61-1;
    int degA=1000;
    int degB=1000;
    vector<LONG> A0 = genVec(degA, p);
    vector<LONG> B0 = genVec(degB, p);

    
    // Addition with O(max(DegA,DegB)) space.:
    benchMark("Addition with extra space: ",2000,[&](){
        auto A=A0;
        auto B=B0;
        auto res=pADDNEW(A,B,degA,degB,p);
    });

    // Addition inplace. This is faster. 
    benchMark("Addition in Place: ",2000,[&](){
        auto A=A0;
        auto B=B0;
        auto res=pADDIP(A,B,degA,degB,p);
    });
    
    // Multiplication:
    /*benchMark("pMul",10,[&](){
        auto res=pMul(A,B,degA,degB,p);
    });*/

    /*
    benchMark("pMulASM2", 10, [&](){
        auto A = A0;           // copy
        auto B = B0;           // copy
        auto res = pMulASM2(A, B, degA, degB, p);
    });


    benchMark("pMulASM",10,[&](){
        auto A = A0;           // copy
        auto B = B0;           // copy
        auto res = pMulASM(A, B, degA, degB, p);
    });

    // Approx. 8000 mu/ops for pMul. 
    // Approx. 13600 mu/ops for pMulASM.

    /*
    // Division:
    benchMark("pDiv",200,[&](){
        vector<LONG> tmp=A;         
        getQuoRemDeg qr=pDiv(tmp,B,degA,degB,p);
        sink ^= qr.degQuo^(qr.degRem<<16);
        if (qr.degRem>=0) sink ^= tmp[qr.degRem];
    });

    // GCD:
    benchMark("pGCD",200,[&](){
        auto g=pGCD(A,B,degA,degB,p);
    });
    */ 
    return 0;
}
