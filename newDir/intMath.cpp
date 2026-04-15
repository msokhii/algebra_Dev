#include<iostream> 
#include<cstdint> 
#include<random> 
#include<vector> 
#include<unordered_map>
#include<algorithm>
#include"int128g.c"

using namespace std;
using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

ULNG seed=1;
ULNG mult=6364136223846793003LL;

// 0<=x<p.

/* 
INTEGER MATH ROUTINES: 
1. RAND64S
2. ADD64S 
3. SUB64S
4. MUL64S
5. MUL64ASM1
6. MUL64ASM2 
7. POWMOD64S
8. INVMOD64S
*/

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

inline LONG add64b(LONG a,LONG b,LONG p){
    LONG r=(a+b)-p;
    r+=(r>>63)&p;
    return r;
}

inline LONG sub64b(LONG a,LONG b,LONG p){
    LONG r=(a-b);
    r+=(r>>63)&p;
    return r;
}

inline LONG mul64bASM(LONG a,LONG b, LONG p){
    ULNG128 res=(ULNG128)a*b;
    ULNG r=(ULNG)(res%p);
    return r;
}

inline LONG neg64s(LONG a,LONG p){ 
    return (a==0)?0:p-a; 
};

// Assuming 0<=a,b<p for the following routines. 

inline LONG mul64bASM(LONG a,LONG b,LONG p){
    LONG q, r;
    __asm__ __volatile__(           \
    "       mulq    %%rdx           \n\t" \
    "       divq    %4              \n\t" \
    : "=a"(q), "=d"(r) : "0"(a), "1"(b), "rm"(p));
    return r;
}

inline LONG mul64bASM2(LONG a,LONG b,LONG p){
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

inline LONG powmod64s(LONG a,LONG n,LONG p){   
    LONG r,s;
    a+=(a>>63)&p; // No bad input.
    if(n==0){return 1;}
    if(n==1){return a;}
    for(r=1,s=a;n>0;n/=2){ 
        if(n&1){
            r=mul64bASM(r,s,p); 
            s=mul64bASM(s,s,p); 
        }
    }
    return r;
};

inline LONG modinv64b(LONG c,LONG p){   
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

/*
POLYNOMIAL ROUTINES: 
*/

struct RatReconFastWS{
    vector<LONG> r1;
    vector<LONG> r2;
    vector<LONG> t1;
    vector<LONG> t2;
    vector<LONG> q;
    vector<LONG> tmpT;

    RatReconFastWS(int degM){
        int n=degM+1;
        r1.resize(n);
        r2.resize(n);
        t1.assign(n,0);
        t2.assign(n,0);
        q.assign(n,0);
        tmpT.assign(n,0);
    }
};

/*
This struct is for pGCDEXTFULL.
*/
struct GCDEX{
	vector<LONG> r;
	vector<LONG> s;
	vector<LONG> t;
	int degR;	
	int degS;
	int degT;
};

/*
This struct is for rational function reconstruction.
*/
struct pairRFR{
	vector<LONG> r;
	vector<LONG> t;
	int degR;
	int degT;
	int flag;
};

/*
This struct returns all values of r,s,t at each iteration.
*/
struct GCDEXHIST{
	GCDEX g;
	vector<vector<LONG>> rTrace;
	vector<vector<LONG>> sTrace;
	vector<vector<LONG>> tTrace;
	vector<int> degRT;
	vector<int> degST;
	vector<int> degTT;
};

/******************************************************************************************/
/* Fast CPU routines                                                                      */
/******************************************************************************************/

#define ZMUL(z,a,b) do { \
    __asm__( \
    "       mulq    %%rdx   \n\t" \
            : "=a"(z[0]), "=d"(z[1]) \
            : "a"(a), "d"(b)); \
} while (0)

#define ZFMA(z,a,b) do { \
    unsigned long u,v; \
    __asm__( \
    "       mulq    %%rdx           \n\t" \
    "       addq    %%rax, %0       \n\t" \
    "       adcq    %%rdx, %1       \n\t" \
            : "=&r"(z[0]), "=&r"(z[1]), "=a"(u), "=d"(v) \
            : "0"(z[0]), "1"(z[1]), "a"(a), "d"(b)); \
} while (0)

#define ZMOD(z,p) __asm__(\
    "       divq    %4      \n\t" \
    "       xorq    %0, %0  \n\t" \
            : "=a"(z[1]), "=d"(z[0]) : "a"(z[0]), "d"(z[1] < p ? z[1] : z[1] % p), "r"(p))

vector<LONG> genVEC64(const int deg,const LONG p){
vector<LONG> v;
/*
Time complexity: O(d+1). 
Space complexity: O(d+1).
Auxillary space: O(1).
*/
if(deg==-1){return v;}
v.resize(deg+1);
for(int i=0;i<=deg;i++){
     v[i]=rand64s(p);
}
return v;
};

vector<LONG> vecCOPY64(const vector<LONG> &v){
vector<LONG> temp; 
temp=v;
return temp;
};

void dispVEC64(const vector<LONG> &v){
if(v.size()==0) cout<<"O"<<"\n";
cout<<"[ ";
for(int i=0;i<v.size();i++){
    if(i==v.size()-1){
        cout<<v[i]<<"*x^"<<i<<"";
        break;
    }
    cout<<v[i]<<"*x^"<<i<<" + "<<"";
    }
    cout<<" ]"<<'\n';
};

// Returns a pair containing the new vector c=(a+b) mod p
// and the degree of c. 

pair<vector<LONG>,int> pADDNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p){
	vector<LONG> c;
	int degC=-1;
	if(degA==-1 && degB==-1) return{c,degC};
	if(degA==-1)return{b,degB};
	if(degB==-1)return{a,degA};
	degC=max(degA,degB);
	c.resize(degC+1,0);
	int i=0;
	while(i<=degA && i<=degB){
		c[i]=add64b(a[i],b[i],p);
		i++;
	}
	for(;i<=degA;i++){
		c[i]=a[i];
	}
	for(;i<=degB;i++){
		c[i]=b[i];
	}
	while(degC>=0 && c[degC]==0) degC--;
	if(degC==-1){
		c.clear();
		return{c,degC};
	}
	c.resize(degC+1);
	return{c,degC};
}

// In place addition. Overwrites a and returns the new degree.

int pADDIP64(vector<LONG> &a,const vector<LONG> &b,int &degA,const int degB,const LONG p){
	if(degA==-1&&degB==-1) return -1;
	if(degB==-1) return degA;
	if(degA==-1){
		a=b;
		degA=degB;
		return degA;
	}
	int maxDeg=max(degA,degB);
	if(a.size()<maxDeg+1) a.resize(maxDeg+1,0);
	int i=0;
	while(i<=degA && i<=degB){
		a[i]=add64b(a[i],b[i],p);
		i++;
	}
	for(;i<=degB;i++) a[i]=b[i];
	while(maxDeg>=0 && a[maxDeg]==0) maxDeg--;
	if(maxDeg==-1){
		a.clear();
		return maxDeg;
	}
	a.resize(maxDeg+1);
	degA=maxDeg;
	return degA;
}

// Returns a pair containing the new vector c=(a-b) mod p
// and the degree of c. 

pair<vector<LONG>,int> pSUBNEW64(const vector<LONG> &a,const vector<LONG> &b,
                                const int degA,const int degB,const LONG p){
    vector<LONG> c;
    int degC=-1;
    if(degA==-1 && degB==-1) return{c,degC};
    if(degA==-1){
        c.resize(degB+1,0);
        for(int i=0;i<=degB;i++) c[i]=neg64s(b[i],p);
        return{c,degB};
    }
    if(degB==-1) return{a,degA};

    degC=max(degA,degB);
    c.resize(degC+1,0);

    int i=0;
    while(i<=degA && i<=degB){
        c[i]=sub64b(a[i],b[i],p);
        i++;
    }
    for(; i<=degA; i++) c[i]=a[i];
    for(; i<=degB; i++) c[i]=neg64s(b[i],p);

    while(degC>=0 && c[degC]==0) degC--;
    if(degC==-1){ c.clear(); return{c,degC}; }
    c.resize(degC+1);
    return{c,degC};
}


// In place subtraction. Overwrites a and returns the new degree.

int pSUBIP64(LONG *a,
             const LONG *b,
             int degA,
             const int degB,
             const LONG p){
	if(degA==-1&&degB==-1) return -1;
	if(degB==-1) return degA;
    if(degA==-1){
		for(int i=0;i<=degB;i++){
			a[i]=neg64s(b[i],p);
		}
		return degB;
	}
	int maxDeg=max(degA,degB);
	int i=0;
	while(i<=degA&&i<=degB){
		a[i]=sub64b(a[i],b[i],p);
		i++;
	}
	for(;i<=degB;i++) a[i]=neg64s(b[i],p);
	while(maxDeg>=0&&a[maxDeg]==0) maxDeg--;
	
	return maxDeg;
}

// Returns a pair containing the new vector c=(a*b) mod p
// and the degree of c. 

pair<vector<LONG>,int> pMULNEW64(const vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p){
	vector<LONG> c;
	if(degA<0 || degB<0) return {c,-1};
	int degC=degA+degB;
	c.resize(degC+1,0);
	for(int i=0;i<=degA;i++){
		for(int j=0;j<=degB;j++){
			LONG prod=mul64bASM(a[i],b[j],p);
			c[i+j]=add64b(c[i+j],prod,p);
		}
	}
	while(degC>=0 && c[degC]==0) degC--;
	if(degC==-1){
		c.clear();
		return {c,degC};
	}
	return {c,degC};
}

// In place multiplication. Overwrites a and returns the new degree.

int pMULIP64(LONG *a,
             const LONG* b,
             int degA,
             int degB,
             const LONG p){
	if(degA<0 || degB<0) return -1;
	int i;
	int k;
	int m;
	int degC=degA+degB; // Called must guarantee that a has enough storage.
	/* 
	If p<2^31 then our product fits inside 2^63 bits.
	We essentially perform a convolution i.e. sum over i of
	a[i]*b[k-i].
	*/
	if(p<2147483648LL){
		LONG t;
		LONG p2; 
		p2=p<<32;
		for(k=degC;k>=0;k--){
			i=max(0,k-degB);
			m=min(k,degA);
			t=0LL; // Accumalator for the running sum.
			while(i<m){
				t-=a[i]*b[k-i];
				i++; 
				t-=a[i]*b[k-i];
				i++;
				t+=(t>>63)&p2;
			}
			if(i==m) t-=a[i]*b[k-i];
			t=(-t)%p;
			t+=(t>>63)&p;
			a[k]=t;
		}
	}
	else{
		ULNG z[2];
		for(k=degC;k>=0;k--){
			i=max(0,k-degB);
			m=min(k,degA);
			z[0]=z[1]=0LL; // 128 bit accumalators.
			while(i<m){
				//
				ZFMA(z,a[i],b[k-i]);
				i++;
				ZFMA(z,a[i],b[k-i]);
				i++;
				if(z[1]>=p) z[1]-=p;
			}
			if(i==m) ZFMA(z,a[i],b[k-i]);
			ZMOD(z,p);
			a[k]=z[0];
		}
	}
	while(degC>=0 && a[degC]==0) degC--;
	return degC;
}

int polMUL64P(LONG *a,
              LONG *b,
              int degA,
              int degB,
              LONG p,
              recint P){
    int i;
    int k;
    int m;
    LONG t;
    if(degA<0 || degB<0){
        return -1;
    }
    int degC=degA+degB;
    for(k=degC;k>=0;k--){
        i=max(0,k-degB);
        m=min(k,degA);
        for(t=0;i<=m;i++){
            t=add64b(t,mulrec64(a[i],b[k-i],P),p);
        }
        a[k]=t;
    }
    while(degC>=0 && a[degC]==0){
        degC--;
    }
    return degC;
}

int polfms64s(LONG *A, LONG *B, LONG *C, int da, int db, int dc, LONG p)
{   // polynomial fused multiply subtract: C -= A*B
    int i,k,m; ULNG z[2];
    if( da<0 || db<0 ) return dc;

    for( k=0; k<=da+db; k++ ) {
        i = max(0, k-db);
        m = min(k, da);
        z[0] = z[1] = 0ll;

        while( i<m ) {
            ZFMA(z, A[i],   B[k-i]); i++;
            ZFMA(z, A[i],   B[k-i]); i++;
            if( z[1] >= p ) z[1] -= p;
        }
        if( i==m ) ZFMA(z, A[i], B[k-i]);

        ZMOD(z, p);

        if( k > dc ) {
            C[k] = (z[0] == 0 ? 0 : p - z[0]);
        } else {
            C[k] = sub64b(C[k], z[0], p);
        }
    }

    for( dc=max(dc, da+db); dc>=0 && C[dc]==0; dc-- );
    return dc;
}



int pMULIP64VANDER(vector<LONG> &a, vector<LONG> &b, int degA, int degB, const LONG p) {
    // Computes: b <- a * b   (in-place on b)
    // a is treated as the left factor (degree degA)
    // b is treated as the right factor and output (degree degB -> degA+degB)

    if (degA < 0 || degB < 0) {
        b.clear();
        return -1;
    }

    if ((int)a.size() < degA + 1) a.resize(degA + 1, 0);
    if ((int)b.size() < degB + 1) b.resize(degB + 1, 0);

    auto normp = [p](LONG x) -> LONG {
        x %= p;
        if (x < 0) x += p;
        return x;
    };

    // Copy inputs (safe even if someone accidentally aliases a and b)
    vector<LONG> A(a.begin(), a.begin() + degA + 1);
    vector<LONG> B(b.begin(), b.begin() + degB + 1);

    for (int i = 0; i <= degA; ++i) A[i] = normp(A[i]);
    for (int j = 0; j <= degB; ++j) B[j] = normp(B[j]);

    int degC = degA + degB;
    vector<LONG> C(degC + 1, 0);

    // Convolution: C[k] = sum_{i=0}^degA A[i] * B[k-i]
    for (int i = 0; i <= degA; ++i) {
        if (A[i] == 0) continue;
        for (int j = 0; j <= degB; ++j) {
            if (B[j] == 0) continue;
            C[i + j] = add64b(C[i + j], mul64bASM(A[i], B[j], p), p);
        }
    }

    // Trim trailing zeros
    while (degC >= 0 && C[degC] == 0) --degC;

    if (degC < 0) {
        b.assign(1, 0);
        return -1;
    }

    C.resize(degC + 1);
    b.swap(C);
    return degC;
}

// Computes scalar polynomial multiplication i.e.
// A=c*A(x) where c is some scalar and returns a new vector.

vector<LONG> polSCMULNEW64(vector<LONG> &a,LONG x,int degA,const LONG p){
	vector<LONG> temp;
	if(x==1){return a;}
	temp.resize(degA+1,0);
	if(x==-1){
		for(int i=0;i<=degA;i++){
			temp[i]=neg64s(a[i],p);
		}
		return temp;
	}
	else{
		for(int i=0;i<=degA;i++){
			temp[i]=mul64bASM(a[i],x,p);
		}
	}
	return temp;	
}

// Computes scalar polynomial multiplication in place i.e.
// A=c*A(x) where c is some scalar. 

void polSCMULIP64(vector<LONG> &a,LONG x,int degA,const LONG p){
	// Since scalar value is 1 no difference.
	if(x==1){return;}
	// If it is negative we just use the neg function.
	if(x==-1){
		for(int i=0;i<=degA;i++){
			a[i]=neg64s(a[i],p);
		}
	}
	else{
		for(int i=0;i<=degA;i++){
			a[i]=mul64bASM(a[i],x,p);
		}
	}
}

// Computes A=A-(ax+b)*B efficiently using accumalators.

int polSUBMUL64(LONG *a,
                const LONG *b,
                LONG aVal,
                LONG bVal,
                int degA,
                int degB,
                const LONG p){
	ULNG z[2];
	LONG t;
	int i;
	/*
	If b is the zero polynomial, then A-(ax+b)*B=A so we simply
	return the degree of A.
	*/
	if(degB<0){
        return degA;
    }
	z[0]=z[1]=0LL;
	/*
	If degA<=degB, then we pad A with zereos. The caller needs to 
    guarantee a[0..degB+1] exists.
	*/
	while(degA<=degB){
    	++degA;
    	a[degA]=0;
    }
	/*
	Constant term is special in the sense b*B does not 
	have any effect on the degrees so we can compute A=b*B directly.
	*/
	t=mul64bASM(bVal,b[0],p);
	a[0]=sub64b(a[0],t,p);
	/*
	Basic for loop using 128 bit accumalators to compute 
	A[i]-(aVal*x-bVal)*B[i].
	*/
	for(i=1;i<=degB;i++){
        z[0]=z[1]=0ULL;
		ZMUL(z,aVal,b[i-1]);
		ZFMA(z,bVal,b[i]);
		ZMOD(z,p);
		t=a[i]-(LONG)z[0];
		a[i]=t+((t>>63)&p);
	}
	/*
	Here, we are updating the new leading coefficient.
	*/
	t=mul64bASM(aVal,b[degB],p);
	a[degB+1]=sub64b(a[degB+1],t,p);
	while(degA>=0 && (a[degA]==0 || a[degA]==p)){
        degA--;
    }
	return degA;
}

int polSUBMUL64P(LONG *a,
                 const LONG *b,
                 LONG aVal,
                 LONG bVal,
                 int degA,
                 int degB,
                 const LONG p,
                 recint P){
    LONG s;
    LONG t;
    int i;
    int d;
    if(degB==-1){
        return degA;
    }
    d=degA;
    while(degA<=degB){
        ++degA;
        a[degA]=0;
    }
    t=mulrec64(bVal,b[0],P);
    a[0]=sub64b(a[0],t,p);
    for(i=1;i<=degB;i++){
        t=mulrec64(aVal,b[i-1],P);
        t=add64b(t,mulrec64(bVal,b[i],P),p);
        a[i]=sub64b(a[i],t,p);
    }
    t=mulrec64(aVal,b[degB],P);
    a[degB+1]=sub64b(a[degB+1],t,p);
    while(degA>=0 && (a[degA]==0 || a[degA]==p)){
        degA--;
    }
    if(degA==d){
        printf("FAIL");
    }
    return degA;
}

// Evaluates a polynomial using Horners rule.

LONG evalHORN64(vector<LONG>& a,LONG alpha,LONG p){
    LONG r = 0LL;
	for (int k=a.size();k-->0;){
        r=add64b(mul64bASM(r,alpha,p),a[k],p);
    }
    return r;
}

LONG pEVAL64(LONG *a,int d,LONG x,const LONG p){
	int i;
	LONG r;
	if(d==-1){return 0;}
	for(r=a[d],i=d-1;i>=0;i--){
		r=add64b(a[i],mul64bASM(x,r,p),p);
	}
	return r;
}

// My implementation of division. Same idea as PDIVIP64 but
// I am also returning the degree of both the quotient 
// and remainder.

pair<int,int> pDIVDEG(vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p){
	/* 
	Changes A in place.
	Suppose degA>=degB>=0. Then, time complexity is O((A-B+1)*(B+1)).
	*/
	if(degB<0){
		cout<<"DIV by 0.\n";
		exit(1);
	}
	if(degA<degB)return {-1,degA};
	if(degB == 0){
    	LONG b0 = b[0] % p; if(b0 < 0) b0 += p;
    	if(b0 == 0){ cout<<"DIV by 0.\n"; exit(1); }
    	LONG inv0 = modinv64b(b0, p);
    	for(int i=0;i<=degA;i++) a[i] = mul64bASM(a[i], inv0, p);
    	return {degA, -1};
	}
	LONG LTB=b[degB];
	LONG invLTB=modinv64b(LTB,p);
	for(int i=degA;i>=degB;i--){
		int k=i-degB;
		LONG LR=a[i];
		if(LR==0){
			a[i]=0; 
			continue;
		}
		LONG prod1=mul64bASM(LR,invLTB,p);
		for(int j=0;j<=degB;j++){
			LONG prod2=mul64bASM(prod1,b[j],p);
			a[k+j]=sub64b(a[k+j],prod2,p);
		}
		a[i]=prod1;
	}
	int degQ=degA-degB;
    while(degQ>0 && a[degB+degQ]==0) --degQ;
	if(degQ==0 && a[degB]==0) degQ=-1;
	int degR=degB-1;
    while(degR>0 && a[degR]==0) --degR;
	if(degR==0 && a[0]==0) degR=-1;
    return {degQ,degR};
}

// We divide a(x) by b(x) and put the remainder in the bottom 
// half of a so a[0...deg(b)-1] and quotient in the top half 
// so a[degb...dega] and return the degree of the remainder.

int polDIVIP64(LONG *a,
               const LONG *b,
               int degA,
               int degB,
               const LONG p){
    int dq;
	int dr;
	int k;
	int j;
	int m; 
	LONG t;
	LONG inv;
    if(degB<0){ 
		cout<<"DIV BY 0.\n"; 
		exit(1); 
	}
    if(degA<degB) return degA;
	/* if (degB == 0) {
    // Divide by constant b0 (must be invertible mod p)
    LONG b0 = b[0] % p; if (b0 < 0) b0 += p;
    if (b0 == 0) { cout << "DIV BY 0.\n"; exit(1); }

    LONG inv0 = modinv64b(b0, p);
    for (int i = 0; i <= degA; i++) a[i] = mul64bASM(a[i], inv0, p);

    // remainder is 0
    return -1;
}*/

	/*
	Special case: If we have degA=degB and we have a monic
	divisor. Since, degrees are the same, the leading term 
	of A is the quotient (call it T). So we do A-(t*B) but 
	the first term is implicitly cancelled so we run a loop from 
	0 to degA-1 and update A accordingly. The degree of the remainder 
	will be degB-1 and we can then chop off the trailing zeroes.
	*/
    if(degA==degB && b[degB]==1){
        t=a[degA];
        for(k=0;k<degA;k++){
            if(b[k]){
				a[k]=sub64b(a[k],mul64bASM(t,b[k],p),p);
			}
		}
		for(dr=degA-1;dr>=0 && a[dr]==0;dr--);
        return dr;
    }
    dq=degA-degB;
    dr=degB-1;
	/* 
	We are working in Zp[x] so inverses are guaranteed to 
	exist. If LC(B) is monic then the inverse is simply 1. Else, 
	we compute the inverse.
	*/
    if(b[degB]==1)inv = 1; 
	else inv=modinv64b(b[degB],p);
	if(p<2147483648LL){ 
	LONG p2;
    p2=p<<32;
	/*
	Same idea as multiplication. If we know p<2^31 we can use 
	accumaltors to make things faster. The indices are 
	extrmely confusing.
	*/
    for(k=degA;k>=0;k--){
        t=a[k];
        m=min(dr,k);
        j=max(0,k-dq);
        while(j<m){
            t-=b[j]*a[k-j+degB]; 
			j++;
            t-=b[j]*a[k-j+degB]; 
			j++;
            t+=(t>>63)&p2;
        }
        if(j==m){
            t-=b[j]*a[k-j+degB];
        }
        t=t%p;
        t+=(t>>63)&p;
        if(k>=degB && inv!=1){
            t=mul64bASM(t,inv,p);
        }
        a[k]=t;
    }
} else{
	ULNG z[2];
    for(k=degA;k>=0;k--){
        z[0]=z[1]=0LL;
        m=min(dr,k);
        j=max(0,k-dq);
        while(j<m){
            ZFMA(z,b[j],a[k-j+degB]); 
			j++;
            ZFMA(z,b[j],a[k-j+degB]); 
			j++;
            if(z[1]>=p)z[1]-=p;
        }
        if(j==m){
            ZFMA(z,b[j],a[k-j+degB]);
        }
        ZMOD(z,p);
        t=a[k]-z[0];
        t+=(t>>63)&p;
        if(k>=degB && inv!=1){
            t=mul64bASM(t,inv,p);
        }
        a[k]=t;
    }
}
    while(dr>=0 && a[dr]==0){
        dr--;
    }
    return dr;
}

int polDIVP(LONG *a,
            LONG *b,
            int degA,
            int degB,
            LONG p,
            recint P){
    int degQ;
    int degR;
    int k;
    int j;
    int m;
    LONG t;
    LONG inv; 
    if(degB<0){
        printf("Div by 0\n");
        return -1;
    }
    if(degA<degB){
        return degA;
    }
    degQ=degA-degB;
    degR=degB-1;
    if(b[degB]==1){
        inv=1;
    }
    else{
        inv=modinv64b(b[degB],p);
    }
    for(k=degA;k>=0;k--){
        t=a[k];
        m=min(degR,k);
        j=max(0,k-degQ);
        for(t=a[k];j<=m;j++){
            t=sub64b(t,mulrec64(b[j],a[k-j+degB],P),p);
        }
        if(k>=degB && inv!=1){
            t=mulrec64(t,inv,P);
        }
        a[k]=t;
    }
    while(degR>0 &&  a[degR]==0){
        degR--;
    }
    return degR;
}

// Makes a polynomial monic. We can do this as we are working 
// in Zp[x] and this is a field so inverses exist.

void polMAKEMONIC64(vector<LONG> &a,const LONG p){
	int degA=a.size()-1;
	if(degA<0 || a[degA]==1) return;
	LONG invTerm;
	invTerm=modinv64b(a[degA],p);
	for(int i=0;i<degA;i++){
		a[i]=mul64bASM(invTerm,a[i],p);
	}
	a[degA]=1;
}

// My version of computing the gcd(a(x),b(x)). This returns 
// the updated vector a with the GCD and its degree.

pair<vector<LONG>,int> polGCDNEW64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p){
	// Division dominates. 
	// Space is O((degA+1)+(degB+1)).
	if(degA==-1){
		polMAKEMONIC64(b,p);
		return {b,degB};
	}
	if(degB==-1){
		polMAKEMONIC64(a,p);
		return {a,degA};
	}
	if(degA<degB){
		swap(a,b);
		swap(degA,degB);
	}
	while(degB!=-1){
		pair<int,int> QR=pDIVDEG(a,b,degA,degB,p);
        int degR=QR.second;
		a.swap(b);       
        degA=degB;
        degB=degR;
    }
	polMAKEMONIC64(a,p);
	return {a,degA};
}

static inline void makeDenMonicIP64(vector<LONG> &num, int degNum,
    vector<LONG> &den, int degDen,
    const LONG p){
if(degDen < 0) return;

LONG lc = den[degDen] % p;
if(lc < 0) lc += p;

if(lc != 1){
LONG inv = modinv64b(lc, p);
polSCMULIP64(den, inv, degDen, p);
if(degNum >= 0){
polSCMULIP64(num, inv, degNum, p);
}
}
}

static inline bool exactDivIP64(vector<LONG> &num, int &degNum,
const vector<LONG> &den, int degDen,
const LONG p){
if(degDen < 0) return false;
if(degNum < degDen) return false;

auto QR = pDIVDEG(num, den, degNum, degDen, p);
int degQ = QR.first;
int degR = QR.second;

if(degR != -1){
return false;
}

if(degQ < 0){
num.clear();
degNum = -1;
return true;
}

// quotient lives in num[degDen .. degDen+degQ]
for(int k = 0; k <= degQ; k++){
num[k] = num[degDen + k];
}
num.resize(degQ + 1);
degNum = degQ;
return true;
}

static inline void makeDenMonicOut64(LONG *num, int degNum,
    LONG *den, int degDen,
    const LONG p,recint P){
if(den[degDen]!=1){
LONG inv=modinv64b(den[degDen],p);
for(int i=0;i<=degDen;i++){
den[i]=mulrec64(den[i],inv,P);
}
for(int i=0;i<=degNum;i++){
num[i]=mulrec64(num[i],inv,P);
}
}
}

int ratReconFastKernelWS(const vector<LONG> &m,
    const vector<LONG> &u,
    int degM,
    int degU,
    int N,
    int D,
    const LONG p,
    RatReconFastWS &W,
    LONG *rOut,
    int &degROut,
    LONG *tOut,
    int &degTOut,
    recint P){

// Copy inputs into workspace
std::copy_n(m.data(), degM + 1, W.r1.data());
std::copy_n(u.data(), degU + 1, W.r2.data());

// Initialize t-sequence:
// r1 = m, r2 = u
// t1 = 0, t2 = 1
W.t2[0] = 1;

int degA  = degM;
int degB  = degU;
int degT1 = -1;
int degT2 = 0;

// Ensure degA >= degB initially
if(degA < degB){
std::swap(W.r1, W.r2);
std::swap(degA, degB);
std::swap(W.t1, W.t2);
std::swap(degT1, degT2);
}

while(degB != -1){

// Stop at the first index k such that deg(r_k) == N
if(degB == N){
    degROut = degB;
    degTOut = degT2;

    std::copy_n(W.r2.data(), degROut + 1, rOut);
    std::copy_n(W.t2.data(), degTOut + 1, tOut);

    // Normalize so denominator is monic
    if(degTOut >= 0){
        LONG lc = tOut[degTOut];
        if(lc == 0){
            degROut = -1;
            degTOut = -1;
            return -30; // unexpected bad denominator
        }

        if(lc != 1){
            LONG lcInv = modinv64b(lc, p);
            if(lcInv == 0){
                degROut = -1;
                degTOut = -1;
                return -31; // inverse does not exist
            }

            for(int i = 0; i <= degROut; i++){
                rOut[i] = mulrec64(rOut[i], lcInv, P);
            }
            for(int i = 0; i <= degTOut; i++){
                tOut[i] = mulrec64(tOut[i], lcInv, P);
            }
        }
    }

    return 0;
}

LONG uInv, aVal, bVal;
int degR, degQ, degT;

// Special degree-1 quotient step
if(degB > 0 && degA - degB == 1){
uInv = modinv64b(W.r2[degB], p);

aVal = mulrec64(W.r1[degA], uInv, P);

bVal = mulrec64(aVal, W.r2[degB - 1], P);
bVal = mulrec64(uInv, sub64b(W.r1[degA - 1], bVal, p), P);

degR = polSUBMUL64P(W.r1.data(), W.r2.data(),
           aVal, bVal, degA, degB, p, P);

degT = polSUBMUL64P(W.t1.data(), W.t2.data(),
           aVal, bVal, degT1, degT2, p, P);
}
else{
// Divide r1 by r2:
// quotient goes into high part of W.r1, remainder stays in low part
degR = polDIVP(W.r1.data(), W.r2.data(), degA, degB, p, P);
degQ = degA - degB;

for(int i = 0; i <= degQ; i++){
W.q[i] = W.r1[degB + i];
}

if(degT2 >= 0){
for(int i = 0; i <= degT2; i++){
W.tmpT[i] = W.t2[i];
}

int degTmpT = degT2;
degTmpT = polMUL64P(W.tmpT.data(), W.q.data(), degTmpT, degQ, p, P);
degT = pSUBIP64(W.t1.data(), W.tmpT.data(), degT1, degTmpT, p);
}
else{
degT = degT1;
}
}

if(degR < 0){
break;
}

// Shift:
// (r1,r2) <- (r2,r)
// (t1,t2) <- (t2,t)
std::swap(W.r1, W.r2);
degA = degB;
degB = degR;

std::swap(W.t1, W.t2);
int oldDegT2 = degT2;
degT2 = degT;
degT1 = oldDegT2;
}

degROut = -1;
degTOut = -1;
return -20;
};

int ratReconNormal(const vector<LONG> &m,
                   const vector<LONG> &u,
                   int degM,
                   int degU,
                   int N,
                   int D,
                   const LONG p,
                   RatReconFastWS &W,
                   LONG *rOut,
                   int &degROut,
                   LONG *tOut,
                   int &degTOut){
    std::copy_n(m.data(),degM+1,W.r1.data());
    std::copy_n(u.data(),degU+1,W.r2.data());
    W.t2[0]=1;
    int degA=degM;
    int degB=degU;
    int degT1=-1;
    int degT2=0;
    if(degA<degB){
        std::swap(W.r1,W.r2);
        std::swap(degA,degB);
        std::swap(W.t1,W.t2);
        std::swap(degT1,degT2);
    }
    while(degB!=-1){
        if(degB==N){
            degROut=degB;
            degTOut=degT2;
            std::copy_n(W.r2.data(),degROut+1,rOut);
            std::copy_n(W.t2.data(),degTOut+1,tOut);
            return 0;
        }
        LONG uInv,aVal,bVal;
        int degR,degQ,degT;
        if(degB>0 && degA-degB==1){
            uInv=modinv64b(W.r2[degB],p);
            aVal=mul64bASM(W.r1[degA],uInv,p);
            bVal=mul64bASM(aVal,W.r2[degB-1],p);
            bVal=mul64bASM(uInv,sub64b(W.r1[degA-1],bVal,p),p);
            degR=polSUBMUL64(W.r1.data(),W.r2.data(),
                            aVal,bVal,degA,degB,p);
            degT=polSUBMUL64(W.t1.data(),W.t2.data(),
                            aVal,bVal,degT1,degT2,p);
        }
        else{
            degR=polDIVIP64(W.r1.data(),W.r2.data(),degA,degB,p);
            degQ=degA-degB;
            for(int i=0;i<=degQ;i++){
                W.q[i]=W.r1[degB+i];
            }
            if(degT2>=0){
                for(int i=0;i<=degT2;i++){
                    W.tmpT[i]=W.t2[i];
                }

            int degTmpT=degT2;
            degTmpT=pMULIP64(W.tmpT.data(),W.q.data(),degTmpT,degQ,p);
            degT=pSUBIP64(W.t1.data(),W.tmpT.data(),degT1,degTmpT,p);
            }
        else{
            degT=degT1;
        }
        }
        if(degR<0){
            break;
        }
        std::swap(W.r1,W.r2);
        degA=degB;
        degB=degR;
        std::swap(W.t1,W.t2);
        int oldDegT2=degT2;
        degT2=degT;
        degT1=oldDegT2;
    }
    degROut=-1;
    degTOut=-1;
    return -20;
};

int ratRecon2(const vector<LONG> &m,
                   const vector<LONG> &u,
                   int degM,
                   int degU,
                   int N,
                   int D,
                   const LONG p,
                   RatReconFastWS &W,
                   LONG *rOut,
                   int &degROut,
                   LONG *tOut,
                   int &degTOut,
                   recint P){
    std::copy_n(m.data(),degM+1,W.r1.data());
    std::copy_n(u.data(),degU+1,W.r2.data());
    W.t2[0]=1;
    int degA=degM;
    int degB=degU;
    int degT1=-1;
    int degT2=0;
    if(degA<degB){
        std::swap(W.r1,W.r2);
        std::swap(degA,degB);
        std::swap(W.t1,W.t2);
        std::swap(degT1,degT2);
    }
    while(degB!=-1){
        if(degB==N){
            degROut=degB;
            degTOut=degT2;
            std::copy_n(W.r2.data(),degROut+1,rOut);
            std::copy_n(W.t2.data(),degTOut+1,tOut);
            return 0;
        }
        LONG uInv,aVal,bVal;
        int degR,degQ,degT;
        if(degB>0 && degA-degB==1){
            uInv=modinv64b(W.r2[degB],p);
            aVal=mulrec64(W.r1[degA],uInv,P);
            bVal=mulrec64(aVal,W.r2[degB-1],P);
            bVal=mulrec64(uInv,sub64b(W.r1[degA-1],bVal,p),P);
            degR=polSUBMUL64(W.r1.data(),W.r2.data(),
                            aVal,bVal,degA,degB,p);
            degT=polSUBMUL64(W.t1.data(),W.t2.data(),
                            aVal,bVal,degT1,degT2,p);
        }
        else{
            degR=polDIVIP64(W.r1.data(),W.r2.data(),degA,degB,p);
            degQ=degA-degB;
            for(int i=0;i<=degQ;i++){
                W.q[i]=W.r1[degB+i];
            }
            if(degT2>=0){
                for(int i=0;i<=degT2;i++){
                    W.tmpT[i]=W.t2[i];
                }

            int degTmpT=degT2;
            degTmpT=pMULIP64(W.tmpT.data(),W.q.data(),degTmpT,degQ,p);
            degT=pSUBIP64(W.t1.data(),W.tmpT.data(),degT1,degTmpT,p);
            }
        else{
            degT=degT1;
        }
        }
        if(degR<0){
            break;
        }
        std::swap(W.r1,W.r2);
        degA=degB;
        degB=degR;
        std::swap(W.t1,W.t2);
        int oldDegT2=degT2;
        degT2=degT;
        degT1=oldDegT2;
    }
    degROut=-1;
    degTOut=-1;
    return -20;
};

int main(){
    LONG p=4294967291;
    LONG a=rand64s(p);
    LONG b=rand64s(p);
    LONG c=rand64s(p);
    cout<<a<<" "<<b<<" "<<c<<"\n";
    LONG f1=add64b(a,b,p);
    LONG f2=sub64b(a,b,p);
    LONG f3M=mul64bASM(a,b,p);
    LONG f3M2=mul64bASM(a,b,p);
    LONG f3M3=mul64bASM2(a,b,p);
    cout<<f1<<" "<<f2<<" "<<f3M<<" "<<f3M2<<" "<<f3M3<<"\n";
    return 0;
}
