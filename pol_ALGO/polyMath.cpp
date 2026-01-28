#include<bits/stdc++.h> 
#include"integerMath.h"
#include"helperF.h"

using namespace std; 

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

int pMul(vector<LONG> &a,vector<LONG> &b,int &degA,int &degB,const LONG p){
	if(degA==-1&&degB==-1){return -1;}
	if(degB==-1){a.clear();a.shrink_to_fit();degA=degB;return degA;}
	if(degA==-1){return -1;}
	int maxDeg=degA+degB;
	vector<LONG> temp(maxDeg+1,0);
	for(int i=0;i<=degA;i++){
		for(int j=0;j<=degB;j++){
			LONG prod=mul64b(a[i],b[j],p);
			temp[i+j]=add64b(temp[i+j],prod,p);
		}
	}
	while(maxDeg>=0 && temp[maxDeg]==0){maxDeg--;}
	if(maxDeg==-1){a.clear();a.shrink_to_fit();degA=maxDeg;return degA;}
	temp.resize(maxDeg+1);
	/*
	vector<T>.swap switches the pointers internally. This is
	O(1). No copying of the elements is done. 
	*/
	//temp.swap(a);
	degA=maxDeg;
	return{degA};
}

int pMulASM(vector<LONG> &a,vector<LONG> &b,int &degA,int &degB,const LONG p){
	if(degA==-1&&degB==-1){return -1;}
	if(degB==-1){a.clear();a.shrink_to_fit();degA=degB;return degA;}
	if(degA==-1){return -1;}
	int maxDeg=degA+degB;
	vector<LONG> temp(maxDeg+1,0);
	for(int i=0;i<=degA;i++){
		for(int j=0;j<=degB;j++){
			LONG prod=mul64bASM(a[i],b[j],p);
			temp[i+j]=add64b(temp[i+j],prod,p);
		}
	}
	while(maxDeg>=0 && temp[maxDeg]==0){maxDeg--;}
	if(maxDeg==-1){a.clear();a.shrink_to_fit();degA=maxDeg;return degA;}
	temp.resize(maxDeg+1);
	/*
	vector<T>.swap switches the pointers internally. This is
	O(1). No copying of the elements is done. 
	*/
	//temp.swap(a);
	degA=maxDeg;
	return{degA};
}

int pMulASM2(vector<LONG> &a,vector<LONG> &b,int &degA,int &degB,const LONG p){
	if(degA==-1&&degB==-1){return -1;}
	if(degB==-1){a.clear();a.shrink_to_fit();degA=degB;return degA;}
	if(degA==-1){return -1;}
	int maxDeg=degA+degB;
	vector<LONG> temp(maxDeg+1,0);
	for(int i=0;i<=degA;i++){
		for(int j=0;j<=degB;j++){
			LONG prod=mul64bASM2(a[i],b[j],p);
			temp[i+j]=add64b(temp[i+j],prod,p);
		}
	}
	while(maxDeg>=0 && temp[maxDeg]==0){maxDeg--;}
	if(maxDeg==-1){a.clear();a.shrink_to_fit();degA=maxDeg;return degA;}
	temp.resize(maxDeg+1);
	/*
	vector<T>.swap switches the pointers internally. This is
	O(1). No copying of the elements is done. 
	*/
	//temp.swap(a);
	degA=maxDeg;
	return{degA};
}

// Horners evaluation.
LONG evalPoly(const vector<LONG> &a,LONG alpha,LONG p){
    LONG r=0;
    for (int i=(int)a.size()-1;i>=0;--i){
        r=add64b(mul64b(r,alpha,p),a[i],p);
    }
    return r;
}

struct getQuoRemDeg{
	int degQuo;
	int degRem;
};

getQuoRemDeg pDiv(vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,LONG p){
	/* 
	Changes A in place.
	Suppose degA>=degB>=0. Then, time complexity is O((A-B+1)*(B+1)).
	*/
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

struct getGCD{
	int degG;
	vector<LONG> g;
};

/* getGCD pGCD(vector<LONG> a,vector<LONG> b,int degA,int degB,const LONG p){
	Division dominates. 
	Space is O((degA+1)+(degB+1)).
	if(degA==-1){
		makeMonic(b,degB,p);
		return {degB,b};
	}
	if(degB==-1){
		makeMonic(a,degA,p);
		return {degA,a};
	}
	if(degA<degB){
		swap(a,b);
		swap(degA,degB);
	}
	while(degB!=-1){
		getQuoRemDeg QR=pDiv(a,b,degA,degB,p);
        int degR=QR.degRem;
		a.swap(b);       
        degA=degB;
        degB=degR;
    }

	// Make monic. 
	if(degA>=0){
        LONG LC=a[degA];
        if(LC!=1){
            LONG invLC=modinv64b(LC,p);
            for(int i=0;i<=degA;i++){
                a[i]=mul64b(a[i],invLC,p);
            }
        }
    }
    return {degA,a};
}
*/

struct fullEx{
	vector<LONG> r;
	vector<LONG> s;
	vector<LONG> t;
	int degR;
	int degS;
	int degT;
};

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
