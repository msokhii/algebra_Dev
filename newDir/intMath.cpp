#include<iostream> 
#include<cstdint> 
#include<random> 
#include<vector> 
#include<unordered_map>
#include<algorithm>
#include<time.h>
#include<chrono>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
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

inline LONG mul64b(LONG a,LONG b, LONG p){
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

/*
NEWTON INTERPOLATION ROUTINES:
*/

int newtonInterpMulRec(LONG* x,
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

int newtonInterpMulNormal(LONG* x,
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
            s=add64b(s,mul64bASM(prod,Y[i],p),p);
	        prod=mul64bASM(prod,sub64b(X[j],X[i],p),p);
        }
        if(prod==0){
            return -1;
        }
        Y[j]=mul64bASM(sub64b(Y[j],s,p),modinv64b(prod,p),p);	
    }
    d=n-1;
    while(d>=0&&y[d]==0){
        d--;
    }
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            Y[j]=sub64b(Y[j],mul64bASM(X[d-i],Y[j+1],p),p);        
	    }
    }
    return d;
}

int mkM(vector<LONG>&m,const vector<LONG> &xs,const LONG p){    
    int degM=0;
    std::vector<LONG>linF(2,0);
    linF[1]=1;
    for(int i=0;i<xs.size();i++){
        linF[0]=(xs[i]==0?0:p-xs[i]);
        degM=pMULIP64(m.data(),linF.data(),degM,1,p);
    };
    return degM;
}

/*
NEWTON INTERPOLATION WRAPPER FOR MAPLE.
*/

extern "C" int cppInterp(int xLen,
    const LONG *xIn,
    int yLen,
    const LONG *yIn,
    const LONG p,
    int outLen,
    LONG *yOut,
    int *degOut)
{
// basic checks
if (!xIn || !yIn || !yOut || !degOut) {
return -1;
}
if (xLen <= 0 || yLen <= 0 || outLen <= 0) {
return -2;
}
if (xLen != yLen) {
return -3;
}
if (outLen < yLen) {
return -4;
}

const int n = xLen;

// initialize outputs
*degOut = -1;
for (int i = 0; i < outLen; ++i) {
yOut[i] = 0;
}

// local working copies since kernel overwrites y
vector<LONG> x(xIn, xIn + n);
vector<LONG> y(yIn, yIn + n);

recint P = recip1(p);

int d = newtonInterpMulRec(x.data(), y.data(), n, p, P);

if (d < 0) {
*degOut = -1;
return d;
}

if (d >= outLen) {
*degOut = -1;
return -5;
}

std::copy_n(y.data(), d + 1, yOut);
*degOut = d;

return 0;
}

/*
RATIONAL RECON. WRAPPER FOR MAPLE.
*/

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
    LONG *nOut,
    int *degNOUT,
    int dOutLen,
    LONG *dOut,
    int *degDOUT)
{

if (!M || !U || !nOut || !dOut || !degNOUT || !degDOUT) return -1;
if (degM < 0 || degU < 0) return -1;
if (mLen <= 0 || uLen <= 0 || nOutLen <= 0 || dOutLen <= 0) return -1;
if (degM >= mLen || degU >= uLen) return -1;

*degNOUT = -1;
*degDOUT = -1;

for (int i = 0; i < nOutLen; ++i) nOut[i] = 0;
for (int i = 0; i < dOutLen; ++i) dOut[i] = 0;

int wsSize = std::max(degM, degU) + 1;

std::vector<LONG> m(M, M + (degM + 1));
std::vector<LONG> u(U, U + (degU + 1));
std::vector<LONG> rTmp(wsSize, 0);
std::vector<LONG> tTmp(wsSize, 0);

int degROut = -1;
int degTOut = -1;
RatReconFastWS W(wsSize);
recint P = recip1(p);

int rc = ratReconFastKernelWS(m,
            u,
            degM,
            degU,
            N,
            D,
            p,
            W,
            rTmp.data(),
            degROut,
            tTmp.data(),
            degTOut,
            P);

if (rc != 0) {
*degNOUT = -1;
*degDOUT = -1;
return rc;
}

if (degROut < 0 || degTOut < 0) return -2;
if (degROut >= wsSize || degTOut >= wsSize) return -3;
if (degROut + 1 > nOutLen) return -4;
if (degTOut + 1 > dOutLen) return -5;
std::copy_n(rTmp.data(), degROut + 1, nOut);
std::copy_n(tTmp.data(), degTOut + 1, dOut);
*degNOUT = degROut;
*degDOUT = degTOut;

return 0;
}

/*
int main(){
    LONG p=4294967291; // This is prevprime(2^32-1) from maple. 
    recint P=recip1(p);
    int degN=5;
    int degD=5;
    const int CALLS=1000; 
    const int ITER=7;

    ofstream logFile("benchMark.txt");
    logFile<<"PRIME -> "<<p<<"\n";
    logFile<<"CALLS -> "<<CALLS<<"\n";
    logFile<<left
        <<setw(10)<<"ITER"
        <<setw(10)<<"degN"
        <<setw(10)<<"degD"
        <<setw(28)<<"avgTimeNewton(mulRec)"
        <<setw(28)<<"avgTimeNewton(mul64)"
        <<setw(28)<<"avgTimeRR(No CPU+mulRec)"
        <<setw(28)<<"avgTimeRR(CPU+mul64)"
        <<setw(28)<<"avgTimeRR(CPU+mulRec)"
        << "\n";

    for(int step=1;step<ITER;step++){
        vector<LONG>n(degN+1,0);
        vector<LONG>d(degD+1,0);
        
        for(int i=0;i<degN+1;i++){
            LONG temp=rand64s(p);
            while(temp==0){
                temp=rand64s(p);
            }
            n[i]=temp;
        }
    
        for(int j=0;j<degD+1;j++){
            LONG temp=rand64s(p);
            while(temp==0){
                temp=rand64s(p);
            }
            d[j]=temp;
        }
        if(d[degD]!=1){
            LONG invTerm;
            invTerm=modinv64b(d[degD],p);
            for(int i=0;i<=degD;i++){
                d[i]=mul64b(invTerm,d[i],p);
            }
            for(int j=0;j<=degN;j++){
                n[j]=mul64b(invTerm,n[j],p);
            }
        }

        vector<LONG>nCopy=n;
        vector<LONG>dCopy=d;
        // We need degN+degD+1 points to interpolate. Here we make the x vector. 
        int m=degN+degD+1;
        vector<LONG>x(m,0);
        for(int i=0;i<m;i++){
            x[i]=i+1;
        }   
        
        vector<LONG>y(m,0);
        for(int i=0;i<m;i++){
            LONG denEval=pEVAL64(d.data(),degD,x[i],p);
            if(denEval==0){
                return -1;
            }
            LONG numEval=pEVAL64(n.data(),degN,x[i],p);
            y[i]=mul64b(numEval,modinv64b(denEval,p),p);
        }

        vector<LONG>yCopy(m,0);
        copy(y.begin(),y.end(),yCopy.begin());
        int degU=newtonInterpMulRec(x.data(),yCopy.data(),m,p,P);
        
        // Timer for newton interpolation using mulrec64 routine.
        auto start=chrono::steady_clock::now();
        for(int i=0;i<CALLS;i++){
            copy(y.begin(),y.end(),yCopy.begin());
            int degU=newtonInterpMulRec(x.data(),yCopy.data(),m,p,P);
        };
        auto stop=chrono::steady_clock::now();
        
        // Timer for newton interpolation using mul64b routine.
        auto newton2Start=chrono::steady_clock::now();
        for(int i=0;i<CALLS;i++){
            copy(y.begin(),y.end(),yCopy.begin());
            int degU=newtonInterpMulNormal(x.data(),yCopy.data(),m,p);
        }
        auto newton2Stop=chrono::steady_clock::now();
        
        // Timer for copying y into y0 for newton interpolation.
        auto cpStart=chrono::steady_clock::now();
        for(int i=0;i<CALLS;i++){
            copy(y.begin(),y.end(),yCopy.begin());
        }
        auto cpStop=chrono::steady_clock::now();
        double cpTotal=chrono::duration<double,std::micro>(cpStop-cpStart).count();
        double total=chrono::duration<double,std::micro>(stop-start).count();
        double newton2Total=chrono::duration<double,std::micro>(newton2Stop-newton2Start).count();
        double avgTimeCp=cpTotal/CALLS;
        double avgTimeNewton=(total/CALLS)-avgTimeCp;
        double avgTimeNewton2=(newton2Total/CALLS)-avgTimeCp;
        vector<LONG>M(m+1,0);
        M[0]=1;
        int degM=mkM(M,x,p);

        RatReconFastWS W(degM);
        RatReconFastWS W2(degM);
        RatReconFastWS W3(degM);
        vector<LONG>rOut(m,0);
        vector<LONG>tOut(m,0);
        vector<LONG>rOut2(m,0);
        vector<LONG>tOut2(m,0);
        vector<LONG>rOut3(m,0);
        vector<LONG>tOut3(m,0);
        int degR=-1;
        int degT=-1;
        int flag=-999;
        int degR2=-1;
        int degT2=-1;
        int flag2=-999;
        int degR3=-1;
        int degT3=-1;
        int flag3=-999;
        int degUCP=degU;
        int degNCP=degN;
        int degDCP=degD;
        int mCP=m;
        int degUCP3=degU;
        int degNCP3=degN;
        int degDCP3=degD;
        int mCP3=m;
        vector<LONG> MCP=M;
        vector<LONG> yCP2=y;
        vector<LONG> MCP3=M;
        vector<LONG> yCP3=y;
        auto start2=chrono::steady_clock::now();
        for(int k=0;k<CALLS;k++){
            flag=ratReconFastKernelWS(M,y,m,degU,
            degN,degD,p,W,rOut.data(),degR,tOut.data(),degT,P);
        }
        auto stop2=chrono::steady_clock::now();
        auto rrNormStart=chrono::steady_clock::now();
        for(int k=0;k<CALLS;k++){
            flag2=ratReconNormal(MCP,yCP2,mCP,degUCP,
            degNCP,degDCP,p,W2,rOut2.data(),degR2,tOut2.data(),degT2);
        }
        auto rrNormStop=chrono::steady_clock::now();
        auto rrNorm2Start=chrono::steady_clock::now();
        for(int k=0;k<CALLS;k++){
            flag3=ratRecon2(MCP3,yCP3,mCP3,degUCP3,
            degNCP3,degDCP3,p,W3,rOut3.data(),degR3,tOut3.data(),degT3,P);
        }
        auto rrNorm2Stop=chrono::steady_clock::now();
        double total2=chrono::duration<double,std::micro>(stop2-start2).count();
        double rrNormTotal=chrono::duration<double,std::micro>(rrNormStop-rrNormStart).count();
        double rrNorm2Total=chrono::duration<double,std::micro>(rrNorm2Stop-rrNorm2Start).count();
        double avgTimeRR=total2/CALLS;
        double avgTimeRRNorm=rrNormTotal/CALLS;
        double avgTimeRRNorm2=rrNorm2Total/CALLS;
        
        logFile<<left<<
                 setw(10)<<step<<
                 setw(10)<<degN<<
                 setw(10)<<degD<<
                 setw(28)<<avgTimeNewton<<
                 setw(28)<<avgTimeNewton2<<
                 setw(28)<<avgTimeRR<<
                 setw(28)<<avgTimeRRNorm<<
                 setw(28)<<avgTimeRRNorm2<<
                 "\n";
                  
        degN*=2;
        degD*=2;
    }
    logFile.close();
    return 0;
}
*/


int main() {
    LONG p = 4294967291;   // prevprime(2^32-1)
    recint P = recip1(p);

    int degN = 5;
    int degD = 5;

    const int CALLS = 1000;
    const int ITER  = 11;

    ofstream logFile("benchMark.txt");
    if (!logFile) {
        cerr << "Could not open benchMark.txt\n";
        return 1;
    }

    logFile << "PRIME -> " << p << "\n";
    logFile << "CALLS -> " << CALLS << "\n";
    logFile << left
            << setw(8)  << "ITER"
            << setw(8)  << "degN"
            << setw(8)  << "degD"
            << setw(24) << "NewtonKernelRec_us"
            << setw(24) << "NewtonKernel64_us"
            << setw(24) << "NewtonWrapCPP_us"
            << setw(24) << "RRKernelFastWS_us"
            << setw(24) << "RRWrapCPP_us"
            << "\n";

    for (int step = 1; step < ITER; ++step) {
        // ------------------------------------------------------------
        // Build random monic rational function n(x)/d(x)
        // ------------------------------------------------------------
        vector<LONG> n(degN + 1, 0);
        vector<LONG> d(degD + 1, 0);

        for (int i = 0; i <= degN; ++i) {
            LONG temp = rand64s(p);
            while (temp == 0) temp = rand64s(p);
            n[i] = temp;
        }

        for (int j = 0; j <= degD; ++j) {
            LONG temp = rand64s(p);
            while (temp == 0) temp = rand64s(p);
            d[j] = temp;
        }

        if (d[degD] != 1) {
            LONG invTerm = modinv64b(d[degD], p);
            for (int i = 0; i <= degD; ++i) d[i] = mul64b(invTerm, d[i], p);
            for (int i = 0; i <= degN; ++i) n[i] = mul64b(invTerm, n[i], p);
        }

        // ------------------------------------------------------------
        // Sample the rational function at m = degN+degD+1 points
        // ------------------------------------------------------------
        int m = degN + degD + 1;

        vector<LONG> x(m, 0);
        for (int i = 0; i < m; ++i) x[i] = i + 1;

        vector<LONG> yVals(m, 0);
        for (int i = 0; i < m; ++i) {
            LONG denEval = pEVAL64(d.data(), degD, x[i], p);
            if (denEval == 0) {
                cerr << "Encountered zero denominator evaluation.\n";
                return 1;
            }
            LONG numEval = pEVAL64(n.data(), degN, x[i], p);
            yVals[i] = mul64b(numEval, modinv64b(denEval, p), p);
        }

        // ------------------------------------------------------------
        // Newton kernel timings
        // We time:
        //   (copy + kernel) - (copy only)
        // so the result is kernel-only.
        // ------------------------------------------------------------
        vector<LONG> yWork(m, 0);

        auto cpStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(yVals.begin(), yVals.end(), yWork.begin());
        }
        auto cpStop = chrono::steady_clock::now();
        double copyOnly_us =
            chrono::duration<double, std::micro>(cpStop - cpStart).count() / CALLS;

        // mulRec kernel
        int degU_rec = -1;
        auto nRecStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(yVals.begin(), yVals.end(), yWork.begin());
            degU_rec = newtonInterpMulRec(x.data(), yWork.data(), m, p, P);
        }
        auto nRecStop = chrono::steady_clock::now();
        double newtonRecWithCopy_us =
            chrono::duration<double, std::micro>(nRecStop - nRecStart).count() / CALLS;
        double newtonKernelRec_us = newtonRecWithCopy_us - copyOnly_us;

        if (degU_rec < 0) {
            cerr << "newtonInterpMulRec failed.\n";
            return 1;
        }

        // Recover coefficient vector Ucoeff from mulRec Newton
        copy(yVals.begin(), yVals.end(), yWork.begin());
        degU_rec = newtonInterpMulRec(x.data(), yWork.data(), m, p, P);
        if (degU_rec < 0) {
            cerr << "newtonInterpMulRec failed while building Ucoeff.\n";
            return 1;
        }
        vector<LONG> Ucoeff(yWork.begin(), yWork.begin() + (degU_rec + 1));

        // mul64 kernel
        int degU_64 = -1;
        auto n64Start = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(yVals.begin(), yVals.end(), yWork.begin());
            degU_64 = newtonInterpMulNormal(x.data(), yWork.data(), m, p);
        }
        auto n64Stop = chrono::steady_clock::now();
        double newton64WithCopy_us =
            chrono::duration<double, std::micro>(n64Stop - n64Start).count() / CALLS;
        double newtonKernel64_us = newton64WithCopy_us - copyOnly_us;

        // ------------------------------------------------------------
        // Newton wrapper timing in C++
        // This is the number you compare to Maple's
        // "Local newton routine timing"
        // ------------------------------------------------------------
        vector<LONG> yOutWrap(m, 0);
        int degOutWrap = -1;

        auto nWrapStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            int rc = cppInterp(
                m, x.data(),
                m, yVals.data(),
                p,
                m, yOutWrap.data(),
                &degOutWrap
            );
            if (rc != 0) {
                cerr << "cppInterp failed with rc = " << rc << "\n";
                return 1;
            }
        }
        auto nWrapStop = chrono::steady_clock::now();
        double newtonWrapCPP_us =
            chrono::duration<double, std::micro>(nWrapStop - nWrapStart).count() / CALLS;

        // ------------------------------------------------------------
        // Build M(x) = prod (x - x_i)
        // ------------------------------------------------------------
        vector<LONG> M(degN + degD + 2, 0);   // size m+1
        M[0] = 1;
        int degM = mkM(M, x, p);              // should be m

        if (degM < 0) {
            cerr << "mkM failed.\n";
            return 1;
        }

        vector<LONG> Mcoeff(M.begin(), M.begin() + (degM + 1));

        // ------------------------------------------------------------
        // RR kernel timing
        // We use Ucoeff (coefficients of interpolant), not yVals.
        // We also subtract copy-only cost to isolate the kernel better.
        // ------------------------------------------------------------
        vector<LONG> Mwork = Mcoeff;
        vector<LONG> Uwork = Ucoeff;

        auto rrCopyStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(Mcoeff.begin(), Mcoeff.end(), Mwork.begin());
            copy(Ucoeff.begin(), Ucoeff.end(), Uwork.begin());
        }
        auto rrCopyStop = chrono::steady_clock::now();
        double rrCopyOnly_us =
            chrono::duration<double, std::micro>(rrCopyStop - rrCopyStart).count() / CALLS;

        RatReconFastWS W(degM);
        vector<LONG> rOut(degN + 1, 0);
        vector<LONG> tOut(degD + 1, 0);

        int flag = -999;
        int degR = -1;
        int degT = -1;

        auto rrKernelStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(Mcoeff.begin(), Mcoeff.end(), Mwork.begin());
            copy(Ucoeff.begin(), Ucoeff.end(), Uwork.begin());

            degR = -1;
            degT = -1;
            flag = ratReconFastKernelWS(
                Mwork,
                Uwork,
                degM,
                degU_rec,
                degN,
                degD,
                p,
                W,
                rOut.data(),
                degR,
                tOut.data(),
                degT,
                P
            );

            if (flag != 0) {
                cerr << "ratReconFastKernelWS failed with rc = " << flag << "\n";
                return 1;
            }
        }
        auto rrKernelStop = chrono::steady_clock::now();
        double rrKernelWithCopy_us =
            chrono::duration<double, std::micro>(rrKernelStop - rrKernelStart).count() / CALLS;
        double rrKernelFastWS_us = rrKernelWithCopy_us - rrCopyOnly_us;

        // ------------------------------------------------------------
        // RR wrapper timing in C++
        // This is the number you compare to Maple's
        // "Local ratrecon routine timing"
        // ------------------------------------------------------------
        vector<LONG> nOutWrap(degN + 1, 0);
        vector<LONG> dOutWrap(degD + 1, 0);
        int degNOutWrap = -1;
        int degDOutWrap = -1;

        auto rrWrapStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            int rc = ratRECON_C(
                degM + 1, degM, Mcoeff.data(),
                degU_rec + 1, degU_rec, Ucoeff.data(),
                degN, degD, p,
                degN + 1, nOutWrap.data(), &degNOutWrap,
                degD + 1, dOutWrap.data(), &degDOutWrap
            );
            if (rc != 0) {
                cerr << "ratRECON_C failed with rc = " << rc << "\n";
                return 1;
            }
        }
        auto rrWrapStop = chrono::steady_clock::now();
        double rrWrapCPP_us =
            chrono::duration<double, std::micro>(rrWrapStop - rrWrapStart).count() / CALLS;

        // ------------------------------------------------------------
        // Log
        // ------------------------------------------------------------
        logFile << left
                << setw(8)  << step
                << setw(8)  << degN
                << setw(8)  << degD
                << setw(24) << newtonKernelRec_us
                << setw(24) << newtonKernel64_us
                << setw(24) << newtonWrapCPP_us
                << setw(24) << rrKernelFastWS_us
                << setw(24) << rrWrapCPP_us
                << "\n";

        degN *= 2;
        degD *= 2;
    }

    logFile.close();
    return 0;
}
