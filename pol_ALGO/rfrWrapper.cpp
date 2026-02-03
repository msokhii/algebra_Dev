#include"polyMath.h"
#include"integerMath.h"
#include<algorithm>
#include<vector>
#include<cstdint> 

using namespace std;
using LONG=int64_t;

extern "C" int ratRECON_C(int mLen,
                         int degM,
                         const LONG *M,
                         int uLen,
                         int degU,
                         const LONG *U,
                         const int N,
                         const int D,
                         const LONG p,
                         int nOutLen,
                         LONG* nOut,
                         int* degNOUT,
                         int dOutLen,
                         LONG* dOut,
                         int *degDOUT){
    if(degM<0 || degU<0){return -1;}
    if(degM>=mLen || degU>=uLen){return -1;}
    vector<LONG> m;
    vector<LONG> u;
    m.assign(M,M+(degM+1));
    u.assign(U,U+(degU+1));
    pairRFR res=ratRecon(m,u,degM,degU,N,D,p);
    if(res.flag!=0){return res.flag;}
    if(res.degR+1>nOutLen){return -1;}
    if(res.degT+1>dOutLen){return -1;}
    for(int i=0;i<nOutLen;i++){
        nOut[i]=0;
    }
    for(int i=0;i<dOutLen;i++){
        dOut[i]=0;
    }
    for(int i=0;i<=res.degR;i++){
        nOut[i]=res.r[i];
    }
    for(int i=0;i<=res.degT;i++){
        dOut[i]=res.t[i];
    }
    *degNOUT=res.degR;
    *degDOUT=res.degT;
    return 0;
}