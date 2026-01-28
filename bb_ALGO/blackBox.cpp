#include<bits/stdc++.h> 
#include"integerMath.h"
#include"polyMath.h"

using namespace std; 

// Black Box for a univariate polynomial in Zp[x].
LONG pBlackBox(const LONG &alpha,const LONG p){
    vector<LONG> f={1,2,3,4,5,6,7,8,9,10};
    LONG r=evalPoly(f,alpha,p);
    return r;
}

LONG getDegBB(const LONG p){
    vector<LONG> m={1};
    vector<LONG> gK;
    int k=0;
    while(true){
        LONG alpha=rand64s(p);
        LONG mEval=evalPoly(m,alpha,p);
        while(mEval==0){
            alpha=rand64s(p);
            mEval=evalPoly(m,alpha,p);
        }
        LONG yK=pBlackBox(alpha,p); 
        LONG gVal=evalPoly(gK,alpha,p); 
        LONG num=sub64b(yK,gVal,p); 
        LONG inv=modinv64b(mEval,p);
        LONG vK=mul64b(num,inv,p); 
        if(vK==0){return (k-1);}
        if(gK.size()<m.size()){gK.resize(m.size(),0);}
        for(int i=0;i<m.size();++i){
            LONG term=mul64b(vK,m[i],p);
            gK[i]=add64b(gK[i],term,p);
        }
        while (!gK.empty()&&gK.back()==0){gK.pop_back();}
        // ---- m := m*(x - alpha) mod p ----
        // In-place update using:
        // r[0] = -alpha*old0
        // r[i] = old_{i-1} - alpha*old_i   for 1<=i<=s-1
        // r[s] = old_{s-1}
        const int s=m.size();
        LONG old0=m[0];
        LONG prev=old0;
        for (int i=1;i<s;++i){
            LONG curr=m[i];
            m[i]=sub64b(prev,mul64b(alpha,curr,p),p);
            prev=curr;
        }
        m.push_back(prev);
        m[0]=sub64b(0,mul64b(alpha,old0,p),p);
        while (!m.empty()&&m.back()==0){m.pop_back();}
        k++;
    }
}
