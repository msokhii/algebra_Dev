#include<vector>
#include<cstdint>
#include<cstdlib>
#include"integerMath.h"
#include"int128g.hpp"

using namespace std;

int newtonInterp(LONG* x,
    LONG* y,
    const int n,
    const LONG p,
    recint P){
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
            s=add64b(s,mulrec64(prod,Y[i],P),p);
	        prod=mulrec64(prod,sub64b(X[j],X[i],p),P);
        }
        if(prod==0){
            return -1;
        }
        Y[j]=mulrec64(sub64b(Y[j],s,p),modinv64b(prod,p),P);	
    }
    d=n-1;
    while(d>=0&&y[d]==0){
        d--;
    }
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            Y[j]=sub64b(Y[j],mulrec64(X[d-i],Y[j+1],P),p);        
	    }
    }
    return d;
}

int newtonInterp2(LONG* x,
    LONG* y,
    const int n,
    const LONG p)
{
if (n < 1) {
return -1;
}

LONG *X = x;
LONG *Y = y;

int d, i, j;
LONG prod, s;
ULONG z[2];   // for ZMUL / ZFMA / ZMOD

for (j = 1; j < n; j++) {
const LONG xj = X[j];
s = Y[0];

// prod = (xj - X[0]) mod p
prod = sub64b(xj, X[0], p);

for (i = 1; i < j; i++) {
// s = (s + prod * Y[i]) mod p
z[0] = (ULONG)s;
z[1] = 0;
ZFMA(z, (ULONG)prod, (ULONG)Y[i]);
ZMOD(z, (ULONG)p);
s = (LONG)z[0];

// prod = (prod * (xj - X[i])) mod p
LONG diff = sub64b(xj, X[i], p);
ZMUL(z, (ULONG)prod, (ULONG)diff);
ZMOD(z, (ULONG)p);
prod = (LONG)z[0];
}

if (prod == 0) {
return -1;
}

// Y[j] = ((Y[j] - s) * modinv64b(prod,p)) mod p
LONG invProd = modinv64b(prod, p);
LONG t = sub64b(Y[j], s, p);
ZMUL(z, (ULONG)t, (ULONG)invProd);
ZMOD(z, (ULONG)p);
Y[j] = (LONG)z[0];
}

d = n - 1;
while (d >= 0 && Y[d] == 0) {
d--;
}

for (i = 1; i <= d; i++) {
const LONG xi = X[d - i];
for (j = d - i; j <= d - 1; j++) {
// Y[j] = (Y[j] - xi * Y[j+1]) mod p
ZMUL(z, (ULONG)xi, (ULONG)Y[j + 1]);
ZMOD(z, (ULONG)p);
Y[j] = sub64b(Y[j], (LONG)z[0], p);
}
}

return d;
}