#include<bits/stdc++.h> 
#include"integerMath.h"
#include"helperF.h"

using namespace std; 

pair<vector<LONG>,int> newtonInterp(vector<LONG> &x,vector<LONG> &y,const int n,const LONG p){
    // if(n<1){return {y.empty(),-1};}
    int d,i,j;
    LONG prod;
    LONG t;
    LONG s;
    vector<LONG> b(y.size(),0); //O(n) extra space?
    b[0]=y[0];
    for(j=1;j<n;j++){
        s=b[0];
        prod=sub64b(x[j],x[0],p);
        for(i=1;i<j;i++){
            s=add64b(s,mul64b(prod,b[i],p),p);
            prod=mul64b(prod,sub64b(x[j],x[i],p),p);
        }
        if(prod==0){exit(1);}
        b[j]=mul64b(sub64b(y[j],s,p),modinv64b(prod,p),p);
    }
    d=n-1;
    while(d>=0&&b[d]==0){d--;}
    /*
    Horners Evaluation. This is done in place so no extra
    memory. 
    */
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            b[j]=sub64b(b[j],mul64b(x[d-i],b[j+1],p),p);
        }
    }
    return {b,d};
}