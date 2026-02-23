#include<bits/stdc++.h>
#include"integerMath.h"
#include"polyMath.h"
#include<vector>
#include<cstdint>

using namespace std;

void vSOLVER64(vector<LONG> &m,vector<LONG> &y,const int n,vector<LONG> &a,vector<LONG> &M,const int shift,const LONG p){
    int i;
    int j;
    LONG u;
    LONG s;
    vector<LONG> A(2,0);

    A[1]=1;
    M[0]=1;
    A[0]=1;

    for(i=0;i<n;i++){
        A[0]=neg64s(m[i],p);
        pMULIP64(A,M,1,i,p);
    }

    for(j=0;j<n;j++){
        A[0]=neg64s(m[j],p);
        i=polDIVIP64(M,A,n,1,p);
        for(i=0;i<n;i++){
            M[i]=M[i+1];
        }
        u=evalHORN64(M,m[j],p);
        if(u==0){
            cout<<"ROOTS ARE NOT DISTINCT.\n";
        }
        u=modinv64b(u,p);
        for(s=0,i=0;i<n;i++){
            s=add64b(s,mul64b(M[i],y[i],p),p);
        }
        s=mul64b(u,s,p);
        if(shift!=0){
            u=modinv64b(m[j],p);
            u=powmod64s(u,shift,p);
            s=mul64b(u,s,p);
        }
        a[j]=s;
        pMULIP64(A,M,1,n-1,p);
    }
    return;
}