#pragma once 

#include<cstdint>
extern long long GLOBALMUL;
using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

LONG rand64s(LONG p);
static inline LONG neg64s(LONG a,LONG p){ 
    return (a==0)?0:p-a; 
};
static inline LONG add64b(LONG a,LONG b,LONG p){
    LONG r=(a+b)-p;
    r+=(r>>63)&p;
    return r;
};
static inline LONG sub64b(LONG a,LONG b,LONG p){
    LONG r=(a-b);
    r+=(r>>63)&p;
    return r;
};
static inline LONG mul64b(LONG a,LONG b,LONG p){
    GLOBALMUL++;
    ULNG128 res=(ULNG128)a*b;
    ULNG r =(ULNG)(res%p);
    return r;
};
LONG mul64bASM(LONG a,LONG b,LONG p);
LONG mul64bASM2(LONG a,LONG b,LONG p);
static inline LONG powmod64s(LONG a,LONG n,LONG p){   
    LONG r,s;
    a+=(a>>63)&p; // No bad input.
    if(n==0){return 1;}
    if(n==1){return a;}
    for(r=1,s=a;n>0;n/=2){ 
        if(n&1){
            r=mul64b(r,s,p); 
            s=mul64b(s,s,p); 
        }
    }
    return r;
};
static inline LONG modinv64b(LONG c,LONG p){   
    LONG d,r,q,r1,c1,d1;
    d=p;
    c1=1;
    d1=0;
    while(d!=0){
        q=c/d;
        r=c-q*d; 
        r1=c1-q*d1;
        c=d;
        c1=d1;
        d=r; 
        d1=r1;
    }
    if(c!=1) return(0);
    if(c1<0) c1+=p;
    return c1;
};
