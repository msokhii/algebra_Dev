#ifndef POLYMATH_H
#define POLYMATH_H

#include"integerMath.h"
#include"helperF.h"

using namespace std;

struct GCDEX{
	vector<LONG> r;
	vector<LONG> s;
	vector<LONG> t;
	int degR;
	int degS;
	int degT;
};

vector<LONG> genVEC64(const int deg,const LONG p);
vector<LONG> vecCOPY64(const vector<LONG> &v);
void dispVEC64(const vector<LONG> &v);
unordered_map<LONG,int> checkPOL64(const vector<LONG> &v,const LONG p);
pair<vector<LONG>,int> pADDNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p);
int pADDIP64(vector<LONG> &a,const vector<LONG> &b,int &degA,const int degB,const LONG p);
pair<vector<LONG>,int> pSUBNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p);
int pSUBIP64(vector<LONG> &a,const vector<LONG> &b,int degA,const int degB,const LONG p);
int pMULIP64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
int polSUBMUL64(vector<LONG> &a,vector<LONG> &b,LONG aVal,LONG bVal,int degA,int degB,const LONG p);
pair<vector<LONG>,int> pMULNEW64(const vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p);
vector<LONG> polSCMULNEW64(vector<LONG> &a,LONG x,int degA,const LONG p);
void polSCMULIP64(vector<LONG> &a,LONG x,int degA,const LONG p);
LONG evalHORN64(const vector<LONG> &a,LONG alpha,const LONG p);
pair<int,int> pDIVDEG(vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p);
int polDIVIP64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
void polMAKEMONIC64(vector<LONG> &a,const LONG p);
pair<vector<LONG>,int> polGCDNEW64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
int polGCD64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
GCDEX pGCDEXFULLSLOW(vector<LONG> &r0,vector<LONG> &r1,int degr0,int degr1,const LONG p);
GCDEX pGCDEXFULLFAST(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);

#endif