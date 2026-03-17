#include<vector>
#include<cstdint>
#include<cstdlib>
#include"integerMath.h"

using namespace std;
// Extra memory b.
/*
int newtonInterp(LONG *x, 
                 LONG *y,
                 const int n,
                 const LONG p,
                 LONG *b){
    int d,i,j;
    LONG prod;
    LONG s;
    b[0]=y[0]; //Newton coeffs. are directly being stored in b.
    for(j=1;j<n;j++){
        s=b[0];
        prod=sub64b(x[j],x[0],p);
        for(i=1;i<j;i++){
            s=add64b(s,mul64b(prod,b[i],p),p);
            prod=mul64b(prod,sub64b(x[j],x[i],p),p);
        }
        if(prod==0){
            exit(1);
        }
        b[j]=mul64b(sub64b(y[j],s,p),modinv64b(prod,p),p);
    }
    /*
    d is the degree of the interpolating polynomial. This 
    can be atmost n-1.
    d=n-1;
    while(d>=0&&b[d]==0){
        d--;
    }
    /*
    Horners Evaluation. This is done in place so no extra
    memory. 
    */
    /*
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            b[j]=sub64b(b[j],mul64b(x[d-i],b[j+1],p),p);
        }
    }
    return d;
}
*/

// No extra memory. Use y as OP.

int newtonInterp(LONG* x,
    LONG* y,
    const int n,
    const LONG p){
    if(n<1){
        return -1;
    }
    LONG *X=x;
    LONG *Y=y;
    int d;
    int i;
    int j;
    LONG prod;
    LONG s;
    for(j=1;j<n;j++){
        const LONG xj=X[j];
        s=Y[0];
        prod=sub64b(xj,X[0],p);
        for(i=1;i<j;i++){
            s=add64b(s,mul64b(prod,Y[i],p),p);
	        prod=mul64b(prod,sub64b(X[j],X[i],p),p);
        }
        if(prod==0){
            return -1;
        }
        Y[j]=mul64b(sub64b(Y[j],s,p),modinv64b(prod,p),p);	
    }
    d=n-1;
    while(d>=0&&y[d]==0){
        d--;
    }
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            Y[j]=sub64b(Y[j],mul64b(X[d-i],Y[j+1],p),p);        
	    }
    }
    return d;
}
