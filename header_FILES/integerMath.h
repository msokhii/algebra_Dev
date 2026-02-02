#ifndef INTEGERMATH_H
#define INTEGERMATH_H

#include<cstdint>

using LONG=int_fast64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

LONG rand64s(LONG p);
LONG neg64s(LONG a,LONG p);
LONG add64b(LONG a,LONG b,LONG p);
LONG sub64b(LONG a,LONG b,LONG p);
LONG mul64b(LONG a,LONG b,LONG p);
LONG mul64bASM(LONG a,LONG b,LONG p);
LONG mul64bASM2(LONG a,LONG b,LONG p);
LONG powmod64s(LONG a,LONG n,LONG p);
LONG modinv64b(LONG c,LONG p);

#endif
