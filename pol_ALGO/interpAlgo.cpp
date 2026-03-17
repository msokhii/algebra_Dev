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
    const LONG p
    ){
    if(n<1){
        return -1;
    }
    int d;
    int i;
    int j;
    LONG prod;
    LONG s;
    for(j=1;j<n;j++){
        s=y[0];
        prod=sub64b(x[j],x[0],p);
        for(i=1;i<j;i++){
            s=add64b(s,mul64b(prod,y[i],p),p);
            
	    prod=mul64b(prod,sub64b(x[j],x[i],p),p);
            
	}
        if(prod==0){
            exit(1);
        }
        y[j]=mul64b(sub64b(y[j],s,p),modinv64b(prod,p),p);
	
}
    d=n-1;
    while(d>=0&&y[d]==0){
        d--;
    }
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            y[j]=sub64b(y[j],mul64b(x[d-i],y[j+1],p),p);
            
	}
    }
    return d;
}


/* int newtonInterp2(LONG* y, const int n, const LONG p) {
    if (n < 1) {
        return -1;
    }

    int d, i, j;
    LONG s, prod;

    // Precompute factorials and inverse factorials mod p.
    std::vector<LONG> fact(n, 0);
    std::vector<LONG> invFact(n, 0);

    fact[0] = 1;
    for (j = 1; j < n; j++) {
        fact[j] = mul64b(fact[j - 1], (LONG)j, p);
    }

    for (j = 0; j < n; j++) {
        invFact[j] = modinv64b(fact[j], p);
        if (invFact[j] == 0) {
            std::exit(1);
        }
    }

    // Step 1: compute Newton coefficients in-place in y.
    for (j = 1; j < n; j++) {
        s = y[0];
        prod = j;   // (x_j - x_0) = j - 0 = j

        for (i = 1; i < j; i++) {
            s = add64b(s, mul64b(prod, y[i], p), p);
            prod = mul64b(prod, (LONG)(j - i), p);   // multiply by (x_j - x_i) = j - i
        }

        y[j] = mul64b(sub64b(y[j], s, p), invFact[j], p);
    }

    // Trim trailing zero coefficients in Newton form.
    d = n - 1;
    while (d >= 0 && y[d] == 0) {
        d--;
    }

    // Step 2: convert Newton form to monomial form in-place.
    for (i = 1; i <= d; i++) {
        LONG xi = (LONG)(d - i);   // because x[d-i] = d-i
        for (j = d - i; j <= d - 1; j++) {
            y[j] = sub64b(y[j], mul64b(xi, y[j + 1], p), p);
        }
    }

    return d;
}

*/
