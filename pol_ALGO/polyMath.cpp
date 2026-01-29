// Compile with g++ -I .Iheader_FILES pol_ALGO/polyMath.cpp pol_ALGO/integerMath.cpp

#include<bits/stdc++.h> 
#include"integerMath.h"
#include"helperF.h"

using namespace std;

/*
This struct is for pGCDEXTNEW.
*/
struct GCDEX{
	vector<LONG> r;
	vector<LONG> s;
	vector<LONG> t;
	int degR;
	int degS;
	int degT;
};

#define ZMUL(z,a,b) __asm__(\
        "       mulq    %%rdx   \n\t" \
                : "=a"(z[0]), "=d"(z[1]) : "a"(a), "d"(b))

#define ZFMA(z,a,b) do {        \
        unsigned long u,v;              \
        __asm__(                        \
        "       mulq    %%rdx           \n\t" \
        "       addq    %%rax, %0       \n\t" \
        "       adcq    %%rdx, %1       \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]), "=a"(u), "=d"(v) : "0"(z[0]), "1"(z[1]), "a"(a), "d"(b));\
        } while (0)

#define ZMOD(z,p) __asm__(\
        "       divq    %4      \n\t" \
        "       xorq    %0, %0  \n\t" \
                : "=a"(z[1]), "=d"(z[0]) : "a"(z[0]), "d"(z[1] < p ? z[1] : z[1] % p), "r"(p))

/******************************************************************************************/
/* Polynomial routines                                                                    */
/******************************************************************************************/

// Generates a vector of size deg+1. Values in the vector
// are generated in random mod p and live in [0..p). 

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
}

// Copies a vector into another vector and then returns it.
// O(degV+1) extra space.
// Use cpp .copy() instead.

vector<LONG> vecCOPY64(const vector<LONG> &v){
	vector<LONG> temp; 
	temp=v;
	return temp;
}

// Prints out the coefficient array of a polynomial from 
// low to high (a0+a1x^1+a2x^2+...+anx^n).

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
}

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

// In place addition. Updates a and returns the new degree.

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

pair<vector<LONG>,int> pSUBNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p){
	vector<LONG> c;
	int degC=-1;
	if(degA==-1 && degB==-1) return{c,degC};
	if(degA==-1){
		c.resize(degB+1,0);
		for(int i=0;i<=degB;i++){
			c[i]=neg64s(b[i],p);
		}
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
	for(;i<=degA;i++){
		c[i]=neg64s(a[i],p);
	}
	for(;i<=degB;i++){
		c[i]=neg64s(b[i],p);
	}
	while(degC>=0 && c[degC]==0) degC--;
	if(degC==-1){
		c.clear();
		return{c,degC};
	}
	c.resize(degC+1);
	return{c,degC};
}

// In place subtraction. Updates a and returns the new degree.

int pSUBIP64(vector<LONG> &a,const vector<LONG> &b,int degA,const int degB,const LONG p){
	if(degA==-1&&degB==-1) return -1;
	if(degB==-1) return degA;
	if(degA==-1){
		a=b;
		degA=degB;
		for(int i=0;i<=degA;i++){
			a[i]=neg64s(b[i],p);
		}
		return degA;
	}
	int maxDeg=max(degA,degB);
	if(a.size()<maxDeg+1) a.resize(maxDeg+1,0);
	int i=0;
	while(i<=degA&&i<=degB){
		a[i]=sub64b(a[i],b[i],p);
		i++;
	}
	for(;i<=degB;i++) a[i]=neg64s(b[i],p);
	while(maxDeg>=0&&a[maxDeg]==0) maxDeg--;
	if(maxDeg==-1){
		a.clear();
		return maxDeg;
	}
	a.resize(maxDeg+1);
	degA=maxDeg;
	return degA;
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
			LONG prod=mul64b(a[i],b[j],p);
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

// In place multiplication. Updates a and returns the new degree.

int pMULIP64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p){
	if(degA<0 || degB<0) return -1;
	int i;
	int k;
	int m;
	int degC=degA+degB;
	if(degC>degA) a.resize(degC+1);
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
	a.resize(degC+1);
	return degC;
}

// Computes A=A-(ax+b)*B efficiently using accumalators.

int polSUBMUL64(vector<LONG> &a,vector<LONG> &b,LONG aVal,LONG bVal,int degA,int degB,const LONG p){
	ULNG z[2];
	LONG t;
	int i;
	/*
	If b is the zero polynomial, then A-(ax+b)*B=A so we simply
	return the degree of A.
	*/
	if(degB==-1){return degA;}
	z[0]=z[1]=0LL;
	/*
	If degA<=degB, then we pad A with zereos.
	*/
	while(degA<=degB){a[++degA]=0;}
	/*
	Constant term is special in the sense b*B does not 
	have any effect on the degrees so we can compute A=b*B directly.
	*/
	t=mul64b(bVal,b[0],p);
	a[0]=sub64b(a[0],t,p);
	/*
	Basic for loop using 128 bit accumalators to compute 
	A[i]-(aVal*x-bVal)*B[i].
	*/
	for(i=1;i<=degB;i++){
		ZMUL(z,aVal,b[i-1]);
		ZFMA(z,bVal,b[i]);
		ZMOD(z,p);
		t=a[i]-z[0];
		a[i]=t+((t>>63)&p);
	}
	/*
	Here, we are updating the new leading coefficient.
	*/
	t=mul64b(aVal,b[degB],p);
	a[degB+1]=sub64b(a[degB+1],t,p);
	// Why do we have an extra checking condition?
	while(degA>=0 && (a[degA]==0 || a[degA]==p)){degA--;}
	return degA;
}

// Evaluates a polynomial using Horners rule.

LONG evalHORN64(const vector<LONG> &a,LONG alpha,const LONG p){
    LONG r=0LL;
    for (int i=a.size()-1;i>=0;i--){
        r=add64b(mul64b(r,alpha,p),a[i],p);
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
	if(degB==0){
		cout<<"DIV by 0.\n";
		exit(1);
	}
	if(degA<degB)return {0,degA};
	LONG LTB=b[degB];
	LONG invLTB=modinv64b(LTB,p);
	for(int i=degA;i>=degB;i--){
		int k=i-degB;
		LONG LR=a[i];
		if(LR==0){
			a[i]=0; 
			continue;
		}
		LONG prod1=mul64b(LR,invLTB,p);
		for(int j=0;j<=degB;j++){
			LONG prod2=mul64b(prod1,b[j],p);
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

int polDIVIP64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p){
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
				a[k]=sub64b(a[k],mul64b(t,b[k],p),p);
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
        if(j==m)t-=b[j]*a[k-j+degB];
        t=t%p;
        t+=(t>>63)&p;
        if(k>=degB && inv!=1)t=mul64b(t,inv,p);
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
        if(j==m)ZFMA(z,b[j],a[k-j+degB]);
        ZMOD(z,p);
        t=a[k]-z[0];
        t+=(t>>63)&p;
        if(k>=degB && inv!=1)t=mul64b(t,inv,p);
        a[k]=t;
    }
}
    while(dr>=0 && a[dr]==0)dr--;
    return(dr);
}

// Makes a polynomial monic. We can do this as we are working 
// in Zp[x] and this is a field so inverses exist.

void polMAKEMONIC64(vector<LONG> &a,const LONG p){
	int degA=a.size();
	if(degA<0 || a[degA]==1) return;
	LONG invTerm;
	invTerm=modinv64b(a[degA],p);
	for(int i=0;i<degA;i++){
		a[i]=mul64b(invTerm,a[i],p);
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

// This computes the gcd(a(x),b(x)) and puts the GCD i.e. g
// in A and returns degree of g. Both original a and b are 
// updated and-or destroyed. 

int polGCD64(vector<LONG> &a, vector<LONG> &b, int degA, int degB, const LONG p) {
    int degR;
    vector<LONG> c;
	vector<LONG> d;
    LONG u;
	LONG aVal;
	LONG bVal;

    // Division by 0 polynomial.
    if (degB<0){
        cout<<"DIV by 0.\n";
        exit(1);
    }

    // Switches pointers internally but destroys a,b.
    c.swap(a);
    d.swap(b);

    // Make sure degC>=degD.
    if(degA<degB){
        swap(c,d);
        swap(degA,degB);
    }

    while(true){
        // Special case: quotient must be linear when degA=degB+1 (and degB>0).
        if(degB>0 && degA-degB==1){
            u=modinv64b(d[degB],p);
            aVal=mul64b(c[degA],u,p);
            bVal=mul64b(aVal,d[degB-1],p);
            bVal=mul64b(u,sub64b(c[degA-1],bVal,p),p); // quotient=ax+b.
			degR=polSUBMUL64(c,d,aVal,bVal,degA,degB,p); // c=c-(ax+b)d.
			if(degR>=degB){cout << "FAIL.\n";}
        }else{
            // General case: compute remainder of c by d (in-place in c)
            degR=polDIVIP64(c,d,degA,degB,p);
        }

        // The remainder is zero -> gcd is d.
        if (degR<0){
            a.swap(d);              // move gcd into a.
            polMAKEMONIC64(a,p);
            return degB;
        }

        // Else continue the algorithm.
        swap(c,d); 
        degA=degB;
        degB=degR;
    }
}


// Make this return monic r. 

/*fullEx pGCDEX(vector<LONG> r0,vector<LONG> r1,int degr0,int degr1,const LONG p){
	vector<LONG> s0{1},s1;
	vector<LONG> t0,t1{1};
	int degS0=0;
	int degS1=-1;
	int degT0=-1;
	int degT1=1;
	while(degr1!=-1){
		int degB=degr1;
		getQuoRemDeg QR=pDiv(r0,r1,degr0,degr1,p);
		int degR=QR.degRem;
		int degQ=QR.degQuo;
		vector<LONG> q;
		if(degQ>=0){
			q.resize(degQ+1);
			for(int k=0;k<=degQ;k++){
				q[k]=r0[degB+k];
			}
		}
		auto [qs1,dqs1]=pMul(q,s1,degQ,degS1,p);
        auto [s2,ds2]=pSub(s0,qs1,degS0,dqs1,p);
        auto [qt1,dqt1]=pMul(q,t1,degQ,degT1,p);
        auto [t2,dt2]=pSub(t0,qt1,degT0,dqt1,p);
		if(degR==-1) r0.clear();
		else r0.resize(degR+1);
		r0.swap(r1);
		degr0=degB;
		degr1=degR;
		s0.swap(s1); degS0 = degS1;
        s1.swap(s2); degS1 = ds2;
        t0.swap(t1); degT0 = degT1;
        t1.swap(t2); degT1 = dt2;
	}
	return{r0,s0,t0,degr0,degS0,degT0};
}
*/


