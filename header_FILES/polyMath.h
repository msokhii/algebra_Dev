#pragma once
#ifndef POLYMATH_H
#define POLYMATH_H

#include"integerMath.h"
#include"helperF.h"
#include<vector>
#include<cstdint> 
#include<unordered_map>

using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

using INT32  = int32_t;
using INT64  = int64_t;
using UINT32 = uint32_t;
using UINT64 = uint64_t;

using namespace std;

struct recint {
  UINT64 s;   /* shift */
  UINT64 v;   /* reciprocal */
  UINT64 d0;  /* divisor shifted up */
  UINT64 d1;
};

struct GCDEX{
	vector<LONG> r;
	vector<LONG> s;
	vector<LONG> t;
	int degR;
	int degS;
	int degT;
};

struct GCDEXHIST{
	GCDEX g;
	vector<vector<LONG>> rTrace;
	vector<vector<LONG>> sTrace;
	vector<vector<LONG>> tTrace;
	vector<int> degRT;
	vector<int> degST;
	vector<int> degTT;
};

struct pairRFR{
	vector<LONG> r;
	vector<LONG> t;
	int degR;
	int degT;
	int flag;
};

vector<LONG> genVEC64(const int deg,const LONG p);
vector<LONG> vecCOPY64(const vector<LONG> &v);
void dispVEC64(const vector<LONG> &v);
unordered_map<LONG,int> checkPOL64(const vector<LONG> &v,const LONG p);
LONG powmodP64(LONG a,LONG n,LONG p,recint P);
pair<vector<LONG>,int> pADDNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p);
int pADDIP64(vector<LONG> &a,const vector<LONG> &b,int &degA,const int degB,const LONG p);
pair<vector<LONG>,int> pSUBNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p);
int pSUBIP64(vector<LONG> &a,const vector<LONG> &b,int degA,const int degB,const LONG p);
int pMULIP64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
int pMULIP64VANDER(vector<LONG> &a, vector<LONG> &b, int degA, int degB, const LONG p);
int polSUBMUL64(vector<LONG> &a,vector<LONG> &b,LONG aVal,LONG bVal,int degA,int degB,const LONG p);
pair<vector<LONG>,int> pMULNEW64(const vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p);
vector<LONG> polSCMULNEW64(vector<LONG> &a,LONG x,int degA,const LONG p);
void polSCMULIP64(vector<LONG> &a,LONG x,int degA,const LONG p);
LONG evalHORN64(vector<LONG> &a,LONG alpha,const LONG p);
LONG pEVAL64(vector<LONG> &a,int d,LONG x,const LONG p);
pair<int,int> pDIVDEG(vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p);
int polDIVIP64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
void polMAKEMONIC64(vector<LONG> &a,const LONG p);
pair<vector<LONG>,int> polGCDNEW64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
int polGCD64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
GCDEX pGCDEXFULLSLOW(vector<LONG> &r0,vector<LONG> &r1,int degr0,int degr1,const LONG p);
GCDEX pGCDEXFULLFAST(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
GCDEXHIST pGCDEXSTORE64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
pairRFR ratRecon(const vector<LONG> &m,const vector<LONG> &u,int degM,int degU,int N,int D,const LONG p);
#endif