#include<iostream>
#include<cstdint>
#include<random>
#include"integerMath.h"

using namespace std;
using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

ULNG seed=1;
ULNG mult=6364136223846793003LL;

// 0<=x<p.

LONG rand64s(LONG p){
    LONG x,y;
    extern ULNG seed,mult;
    seed=mult*seed;
    x=seed>>32;
    seed=mult*seed;
    y=seed>>32;
    x=(x<<31) | y;
    x=x%p;
    return(x);
}

// Assuming 0<=a,b<p for the following routines. 

LONG add64b(LONG a,LONG b,LONG p){
    LONG r=(a+b)-p;
    r+=(r>>63)&p;
    return r;
}

LONG sub64b(LONG a,LONG b,LONG p){
    LONG r=(a-b);
    r+=(r>>63)&p;
    return r;
}
 
LONG mul64b(LONG a,LONG b,LONG p){
    ULNG128 res=(ULNG128)a*b;
    ULNG r =(ULNG)(res%p);
    return r;
}    

LONG mul64bASM(LONG a,LONG b,LONG p){
    LONG q, r;
    __asm__ __volatile__(           \
    "       mulq    %%rdx           \n\t" \
    "       divq    %4              \n\t" \
    : "=a"(q), "=d"(r) : "0"(a), "1"(b), "rm"(p));
    return r;
}

LONG mul64bASM2(LONG a,LONG b,LONG p){
    LONG q;
    LONG r;
    __asm__ __volatile__(
        "movq %[p],%%r8\n\t"
        "movq %[a],%%rax\n\t"
        "mulq %[b]\n\t"
        "divq %%r8\n\t"
        :"=&a"(q),"=&d"(r)
        :[a] "r"(a),[b] "rm"(b),[p] "r"(p)
        :"r8","cc"
    );
    return r;
}

LONG neg64s(LONG a,LONG p){ 
    return (a==0)?0:p-a; 
}

LONG powmod64s(LONG a,LONG n,LONG p){   
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
}

LONG modinv64b(LONG c,LONG p){   
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
}


