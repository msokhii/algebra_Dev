#pragma once
#ifndef FASTVSOLVE_H
#define FASTVSOLVE_H

#include<cstdint>
#include<vector>
#include"integerMath.h"
#include"polyMath.h"

using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

void vSOLVER64(vector<LONG> &m,vector<LONG> &y,const int n,vector<LONG> &a,vector<LONG> &M,const int shift,const LONG p);

#endif