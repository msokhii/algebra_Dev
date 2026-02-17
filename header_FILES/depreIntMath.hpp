#pragma once
/*
  int128_arith.hpp  (header-only, C++17)

  - Multiword add/sub/mul/div (128/256-bit) with optional x86-64 GCC inline asm
  - Reciprocal-division routines (Möller & Granlund style)
  - Optional __int128 extensions

  Notes:
  - The 256/128 division path (_div42) assumes little-endian layout when casting
    UINT64* to UINT32*. This matches typical x86-64.
  - Requires C++17 for best behavior, but most is fine earlier too.
*/

#include <cstdint>
#include <cstdlib>
#include <cstring>

#ifndef INT128_ASSEMBLY
#define INT128_ASSEMBLY 1   /* use x86-64 assembly code */
#endif

#ifndef INT128_EXTENDED
#define INT128_EXTENDED 1   /* functions for __int128_t */
#endif

#ifndef INT128_MACROS
#define INT128_MACROS 0     /* keep internal macros defined after include */
#endif

// ----- types -----
using INT32  = std::int32_t;
using INT64  = std::int64_t;
using UINT32 = std::uint32_t;
using UINT64 = std::uint64_t;

static_assert(sizeof(UINT64) == 8, "UINT64 must be 64-bit.");
static_assert(sizeof(UINT32) == 4, "UINT32 must be 32-bit.");

// ----- divide exception -----
[[noreturn]] inline void int128_div_exception() {
#if defined(__GNUC__) || defined(__clang__)
  __builtin_trap();
#else
  std::abort();
#endif
}
#define DIV_EXCEPTION() ::int128_div_exception()

// ----- API -----
void add128(UINT64 *x, const UINT64 *y);                 /* x += y (128-bit) */
void add256(UINT64 *x, const UINT64 *y);                 /* x += y (256-bit) */
void sub128(UINT64 *x, const UINT64 *y);                 /* x -= y (128-bit) */
void sub256(UINT64 *x, const UINT64 *y);                 /* x -= y (256-bit) */

void mul128(UINT64 *x, UINT64 a, UINT64 b);              /* x = a*b (128-bit) */
void mul256(UINT64 *x, const UINT64 *y, const UINT64 *z);/* x = y*z (256-bit) */

void div128(UINT64 *x, UINT64 a);                        /* x = [x/a : x%a] */
void div256(UINT64 *x, const UINT64 *y);                 /* x = [x/y : x%y] */

UINT64 mulmod64(UINT64 a, UINT64 b, UINT64 p);
void   mulmod128(UINT64 *x, const UINT64 *y, const UINT64 *z, const UINT64 *p);

struct recint {
  UINT64 s;   /* shift */
  UINT64 v;   /* reciprocal */
  UINT64 d0;  /* divisor shifted up */
  UINT64 d1;
};

recint recip1(UINT64 p);           /* reciprocal of  64-bit p */
recint recip2(const UINT64 *p);    /* reciprocal of 128-bit p */

void   rec128(UINT64 *x, recint v);                          /* reciprocal div128 */
void   rec256(UINT64 *x, recint v);                          /* reciprocal div256 */

UINT64 mulrec64(UINT64 a, UINT64 b, recint v);
void   mulrec128(UINT64 *x, const UINT64 *y, const UINT64 *z, recint v);

inline void fma128(UINT64 *x, UINT64 a, UINT64 b);           /* x += a*b (128) */
inline void fms128(UINT64 *x, UINT64 a, UINT64 b);           /* x -= a*b (128) */
inline void fma256(UINT64 *x, const UINT64 *y, const UINT64 *z); /* x += y*z (256) */
inline void fms256(UINT64 *x, const UINT64 *y, const UINT64 *z); /* x -= y*z (256) */

#if INT128_EXTENDED
  using INT128  = __int128_t;
  using UINT128 = __uint128_t;
  void    mul256i(UINT64 *x, UINT128 a, UINT128 b);
  void    div256i(UINT64 *x, UINT128 a, UINT128 *q, UINT128 *r);
  UINT128 mulmod128i(UINT128 a, UINT128 b, UINT128 p);
  UINT128 mulrec128i(UINT128 a, UINT128 b, recint v);
#endif

/* -------------------- low-level macros -------------------- */

#if INT128_ASSEMBLY && defined(__GNUC__) && defined(__x86_64__)

  #define __ADD(x,y) "  addq  " #x ", " #y "\n\t"
  #define __ADC(x,y) "  adcq  " #x ", " #y "\n\t"
  #define __SUB(x,y) "  subq  " #x ", " #y "\n\t"
  #define __SBB(x,y) "  sbbq  " #x ", " #y "\n\t"
  #define __MUL(x)   "  mulq  " #x "\n\t"
  #define __DIV(x)   "  divq  " #x "\n\t"

  #define ADD22(x0,x1,y0,y1) __asm__( \
      __ADD(%2, %0) \
      __ADC(%3, %1) \
    : "=r"(x0), "=r"(x1) : "g"(y0), "g"(y1), "0"(x0), "1"(x1) : "cc")

  #define ADD33(x0,x1,x2,y0,y1,y2) __asm__( \
      __ADD(%3, %0) \
      __ADC(%4, %1) \
      __ADC(%5, %2) \
    : "=r"(x0), "=r"(x1), "=r"(x2) \
    : "g"(y0), "g"(y1), "g"(y2), "0"(x0), "1"(x1), "2"(x2) : "cc")

  #define ADD44(x0,x1,x2,x3,y0,y1,y2,y3) __asm__( \
      __ADD(%4, %0) \
      __ADC(%5, %1) \
      __ADC(%6, %2) \
      __ADC(%7, %3) \
    : "=r"(x0), "=r"(x1), "=r"(x2), "=r"(x3) \
    : "g"(y0), "g"(y1), "g"(y2), "g"(y3), \
      "0"(x0), "1"(x1), "2"(x2), "3"(x3) : "cc")

  #define SUB22(x0,x1,y0,y1) __asm__( \
      __SUB(%2, %0) \
      __SBB(%3, %1) \
    : "=r"(x0), "=r"(x1) : "g"(y0), "g"(y1), "0"(x0), "1"(x1) : "cc")

  #define SUB33(x0,x1,x2,y0,y1,y2) __asm__( \
      __SUB(%3, %0) \
      __SBB(%4, %1) \
      __SBB(%5, %2) \
    : "=r"(x0), "=r"(x1), "=r"(x2) \
    : "g"(y0), "g"(y1), "g"(y2), "0"(x0), "1"(x1), "2"(x2) : "cc")

  #define SUB44(x0,x1,x2,x3,y0,y1,y2,y3) __asm__( \
      __SUB(%4, %0) \
      __SBB(%5, %1) \
      __SBB(%6, %2) \
      __SBB(%7, %3) \
    : "=r"(x0), "=r"(x1), "=r"(x2), "=r"(x3) \
    : "g"(y0), "g"(y1), "g"(y2), "g"(y3), \
      "0"(x0), "1"(x1), "2"(x2), "3"(x3) : "cc")

  #define MUL211(x0,x1,a,b) __asm__( \
      __MUL(%3) \
    : "=a"(x0), "=d"(x1) : "0"(a), "r"(b) : "cc")

  #define DIV21H(x0,x1,a) __asm__( \
      __DIV(%2) \
    : "=a"(x0), "=d"(x1) : "r"(a), "0"(x0), "1"(x1) : "cc")

#else
  // generic C/C++
  #define ADD22(x0,x1,y0,y1) do { \
    (x0) += (y0); \
    (x1) += (y1) + ((x0) < (y0)); \
  } while (0)

  #define ADD33(x0,x1,x2,y0,y1,y2) do { \
    UINT64 s1; \
    s1 = (x1); \
    (x0) += (y0); \
    (x1) += ((x0) < (y0)); \
    (x2) += ((x1) < s1); \
    (x1) += (y1); \
    (x2) += (y2) + ((x1) < (y1)); \
  } while (0)

  #define ADD44(x0,x1,x2,x3,y0,y1,y2,y3) do { \
    UINT64 s1,s2; \
    s1 = (x1); s2 = (x2); \
    (x0) += (y0); \
    (x1) += ((x0) < (y0)); \
    (x2) += ((x1) < s1); \
    (x3) += ((x2) < s2); \
    s2 = (x2); \
    (x1) += (y1); \
    (x2) += ((x1) < (y1)); \
    (x3) += ((x2) < s2); \
    (x2) += (y2); \
    (x3) += (y3) + ((x2) < (y2)); \
  } while (0)

  #define SUB22(x0,x1,y0,y1) do { \
    UINT64 s0; \
    s0 = (x0); \
    (x0) -= (y0); \
    (x1) -= (y1) + (s0 < (y0)); \
  } while (0)

  #define SUB33(x0,x1,x2,y0,y1,y2) do { \
    UINT64 s0,s1; \
    s0 = (x0); s1 = (x1); \
    (x0) -= (y0); \
    (x1) -= (s0 < (y0)); \
    (x2) -= (s1 < (x1)); \
    s1 = (x1); \
    (x1) -= (y1); \
    (x2) -= (y2) + (s1 < (y1)); \
  } while (0)

  #define SUB44(x0,x1,x2,x3,y0,y1,y2,y3) do { \
    UINT64 s0,s1,s2; \
    s0 = (x0); s1 = (x1); s2 = (x2); \
    (x0) -= (y0); \
    (x1) -= (s0 < (y0)); \
    (x2) -= (s1 < (x1)); \
    (x3) -= (s2 < (x2)); \
    s1 = (x1); s2 = (x2); \
    (x1) -= (y1); \
    (x2) -= (s1 < (y1)); \
    (x3) -= (s2 < (x2)); \
    s2 = (x2); \
    (x2) -= (y2); \
    (x3) -= (y3) + (s2 < (y2)); \
  } while (0)

  #define MUL211(x0,x1,a,b) do { \
    UINT64 a0,a1,b0,b1,t0,t1; \
    a0 = (a) & 0xFFFFFFFFu; \
    b0 = (b) & 0xFFFFFFFFu; \
    (x0) = a0*b0; \
    a1 = (a) >> 32; \
    b1 = (b) >> 32; \
    (x1) = a1*b1; \
    t0 = a0*b1; \
    t1 = t0 >> 32; \
    t0 <<= 32; \
    ADD22((x0),(x1),t0,t1); \
    t0 = a1*b0; \
    t1 = t0 >> 32; \
    t0 <<= 32; \
    ADD22((x0),(x1),t0,t1); \
  } while (0)
#endif

/* ----------------- multiply-add helpers ----------------- */

#define FMA211(x0,x1,a,b) do { \
  UINT64 u0,u1; \
  MUL211(u0,u1,(a),(b)); \
  ADD22((x0),(x1),u0,u1); \
} while (0)

#define FMS211(x0,x1,a,b) do { \
  UINT64 u0,u1; \
  MUL211(u0,u1,(a),(b)); \
  SUB22((x0),(x1),u0,u1); \
} while (0)

#define MUL321(x0,x1,x2,y0,y1,a) do { \
  (x2) = 0; \
  MUL211((x0),(x1),(y0),(a)); \
  FMA211((x1),(x2),(y1),(a)); \
} while (0)

#define FMA321(x0,x1,x2,y0,y1,a) do { \
  UINT64 u2,u3,u4; \
  MUL321(u2,u3,u4,(y0),(y1),(a)); \
  ADD33((x0),(x1),(x2),u2,u3,u4); \
} while (0)

#define FMS321(x0,x1,x2,y0,y1,a) do { \
  UINT64 u2,u3,u4; \
  MUL321(u2,u3,u4,(y0),(y1),(a)); \
  SUB33((x0),(x1),(x2),u2,u3,u4); \
} while (0)

#define MUL422(x0,x1,x2,x3,y0,y1,z0,z1) do { \
  UINT64 t2,t3; \
  MUL211((x0),(x1),(y0),(z0)); \
  MUL211((x2),(x3),(y1),(z1)); \
  MUL211(t2,t3,(y1),(z0)); \
  ADD33((x1),(x2),(x3),t2,t3,0); \
  MUL211(t2,t3,(y0),(z1)); \
  ADD33((x1),(x2),(x3),t2,t3,0); \
} while (0)

#define FMA422(x0,x1,x2,x3,y0,y1,z0,z1) do { \
  UINT64 u5,u6,u7,u8; \
  MUL422(u5,u6,u7,u8,(y0),(y1),(z0),(z1)); \
  ADD44((x0),(x1),(x2),(x3),u5,u6,u7,u8); \
} while (0)

#define FMS422(x0,x1,x2,x3,y0,y1,z0,z1) do { \
  UINT64 u5,u6,u7,u8; \
  MUL422(u5,u6,u7,u8,(y0),(y1),(z0),(z1)); \
  SUB44((x0),(x1),(x2),(x3),u5,u6,u7,u8); \
} while (0)

/* ------------------- core ops ------------------- */

inline void add128(UINT64 *x, const UINT64 *y) { ADD22(x[0],x[1],y[0],y[1]); }
inline void sub128(UINT64 *x, const UINT64 *y) { SUB22(x[0],x[1],y[0],y[1]); }
inline void mul128(UINT64 *x, UINT64 a, UINT64 b) { MUL211(x[0],x[1],a,b); }

inline void fma128(UINT64 *x, UINT64 a, UINT64 b) { FMA211(x[0],x[1],a,b); }
inline void fms128(UINT64 *x, UINT64 a, UINT64 b) { FMS211(x[0],x[1],a,b); }

inline void add256(UINT64 *x, const UINT64 *y) { ADD44(x[0],x[1],x[2],x[3],y[0],y[1],y[2],y[3]); }
inline void sub256(UINT64 *x, const UINT64 *y) { SUB44(x[0],x[1],x[2],x[3],y[0],y[1],y[2],y[3]); }
inline void mul256(UINT64 *x, const UINT64 *y, const UINT64 *z) { MUL422(x[0],x[1],x[2],x[3],y[0],y[1],z[0],z[1]); }

inline void fma256(UINT64 *x, const UINT64 *y, const UINT64 *z) { FMA422(x[0],x[1],x[2],x[3],y[0],y[1],z[0],z[1]); }
inline void fms256(UINT64 *x, const UINT64 *y, const UINT64 *z) { FMS422(x[0],x[1],x[2],x[3],y[0],y[1],z[0],z[1]); }

/* ------------------- division ------------------- */

/* 128/64 division from Hacker's Delight (fixed for s==0 to avoid >>64 UB) */
static inline void _div21(UINT64 u0, UINT64 u1, UINT64 v, UINT64 *q, UINT64 *r)
{
#ifdef DIV21H
  DIV21H(u0,u1,v);
  *q = u0;
  *r = u1;
#else
  const UINT64 b = (UINT64)1 << 32;

  if (u1 > v) DIV_EXCEPTION(); /* quotient overflow */

  INT64 s = 0;
  while ((v >> 63) == 0) { v <<= 1; ++s; }

  const UINT64 vn1 = v >> 32;
  const UINT64 vn0 = v & 0xFFFFFFFFu;

  const UINT64 un32 = (u1 << s) | (s ? (u0 >> (64 - s)) : 0);
  const UINT64 un10 = u0 << s;
  const UINT64 un1  = un10 >> 32;
  const UINT64 un0  = un10 & 0xFFFFFFFFu;

  UINT64 q1 = un32 / vn1;
  UINT64 rhat = un32 - q1 * vn1;

again1:
  if (q1 >= b || q1*vn0 > b*rhat + un1) {
    q1--;
    rhat += vn1;
    if (rhat < b) goto again1;
  }

  const UINT64 un21 = un32*b + un1 - q1*v;
  UINT64 q0 = un21 / vn1;
  rhat = un21 - q0 * vn1;

again2:
  if (q0 >= b || q0*vn0 > b*rhat + un0) {
    q0--;
    rhat += vn1;
    if (rhat < b) goto again2;
  }

  *r = (un21*b + un0 - q0*v) >> s;
  *q = q1*b + q0;
#endif
}

inline void div128(UINT64 *x, UINT64 a) { _div21(x[0],x[1],a,&x[0],&x[1]); }

inline UINT64 mulmod64(UINT64 a, UINT64 b, UINT64 p)
{
  UINT64 x0,x1;
  MUL211(x0,x1,a,b);
  _div21(x0,x1,p,&x0,&x1);
  return x1;
}

/* multiword division from Hacker's Delight
   this header version is stack-fixed for the only needed sizes: m<=8, n<=4 */
static inline void _divmnu(UINT32 *u, UINT32 *v, UINT32 *q, UINT32 *r, INT32 m, INT32 n)
{
  constexpr INT32 MMAX = 8;
  constexpr INT32 NMAX = 4;
  if (m > MMAX || n > NMAX) DIV_EXCEPTION();

  const UINT64 b = 4294967296ULL;

  UINT32 un[MMAX + 1] = {0};
  UINT32 vn[NMAX]     = {0};

  UINT64 qhat, rhat, p;
  INT64 t, k;
  INT32 s, i, j;

  while (m && u[m-1]==0) m--;
  while (n && v[n-1]==0) n--;
  if (n==0) DIV_EXCEPTION();

  if (n==1) {
    for (k=0, j=m-1; j >= 0; j--) {
      q[j] = (UINT32)((k*b + u[j]) / v[0]);
      k    = (k*b + u[j]) - (UINT64)q[j]*v[0];
    }
    r[0] = (UINT32)k;
    return;
  }

  for (s=0, t=v[n-1]; !(t >> 31); t <<= 1, s++);

  for (i=n-1; i > 0; i--) {
    vn[i] = (UINT32)((v[i] << s) | ((UINT64)v[i-1] >> (32-s)));
  }
  vn[0] = (UINT32)(v[0] << s);

  un[m] = (UINT32)((UINT64)u[m-1] >> (32-s));
  for (i=m-1; i > 0; i--) {
    un[i] = (UINT32)((u[i] << s) | ((UINT64)u[i-1] >> (32-s)));
  }
  un[0] = (UINT32)(u[0] << s);

  for (j=m-n; j >= 0; j--) {
    p = (UINT64)un[j+n]*b + un[j+n-1];
    qhat = p / vn[n-1];
    rhat = p % vn[n-1];
again:
    if (qhat >= b || qhat*vn[n-2] > b*rhat + un[j+n-2]) {
      qhat--;
      rhat += vn[n-1];
      if (rhat < b) goto again;
    }
    for (k=0, i=0; i < n; i++) {
      p = qhat*vn[i];
      t = (INT64)un[i+j] - (INT64)k - (INT64)(p & 0xFFFFFFFFULL);
      un[i+j] = (UINT32)t;
      k = (p >> 32) - ((UINT64)t >> 32);
    }
    t = (INT64)un[j+n] - (INT64)k;
    un[j+n] = (UINT32)t;
    q[j] = (UINT32)qhat;
    if (t < 0) {
      q[j] = (UINT32)(q[j] - 1);
      for (k=0, i=0; i < n; i++) {
        t = (UINT64)un[i+j] + vn[i] + (UINT64)k;
        un[i+j] = (UINT32)t;
        k = t >> 32;
      }
      un[j+n] = (UINT32)(un[j+n] + k);
    }
  }

  for (i=0; i < n-1; i++) {
    r[i] = (UINT32)((un[i] >> s) | ((UINT64)un[i+1] << (32-s)));
  }
  r[n-1] = (UINT32)(un[n-1] >> s);
}

/* 256/128 long division */
static inline void _div42(UINT64 *x, const UINT64 *y, UINT64 *q, UINT64 *r)
{
  UINT32 R[10] = {0};
  UINT32 Q[10] = {0};

  if (x[3] > y[1] || (x[3]==0 && y[1]==0 && x[2] > y[0])) DIV_EXCEPTION();

  if (x[3]==0 && x[2]==0 && (x[1] < y[1] || (x[1]==y[1] && x[0] < y[0]))) {
    q[0]=q[1]=0;
    r[0]=x[0]; r[1]=x[1];
    return;
  }

  _divmnu(reinterpret_cast<UINT32*>(x),
          const_cast<UINT32*>(reinterpret_cast<const UINT32*>(y)),
          Q, R, 8, 4);

  r[0] = (UINT64(R[1]) << 32) | R[0];
  r[1] = (UINT64(R[3]) << 32) | R[2];
  q[0] = (UINT64(Q[1]) << 32) | Q[0];
  q[1] = (UINT64(Q[3]) << 32) | Q[2];
}

inline void div256(UINT64 *x, const UINT64 *y)
{
  UINT64 q[2], r[2];
  _div42(x, y, q, r);
  x[0]=q[0]; x[1]=q[1];
  x[2]=r[0]; x[3]=r[1];
}

inline void mulmod128(UINT64 *x, const UINT64 *y, const UINT64 *z, const UINT64 *p)
{
  UINT64 t[4], q[2], r[2];
  MUL422(t[0],t[1],t[2],t[3],y[0],y[1],z[0],z[1]);
  _div42(t, p, q, r);
  x[0]=r[0]; x[1]=r[1];
}

/* ---------------- inverse division (reciprocal) ---------------- */

static inline UINT64 _reciprocal21(UINT64 p)
{
  UINT64 u0 = ~UINT64(0);
  UINT64 u1 = -(p+1);
  UINT64 q, r;
  _div21(u0,u1,p,&q,&r);
  return q;
}

static inline UINT64 _reciprocal32(UINT64 d0, UINT64 d1)
{
  UINT64 t0, t1, p, v;
  v = _reciprocal21(d1);
  p = d1*v + d0;
  if (p < d0) {
    v--;
    if (p >= d1) { v--; p -= d1; }
    p -= d1;
  }
  MUL211(t0,t1,v,d0);
  p += t1;
  if (p < t1) {
    v--;
    if (p >= d1 || (p==d1 && t0 >= d0)) v--;
  }
  return v;
}

inline recint recip1(UINT64 p)
{
  if (p==0) DIV_EXCEPTION();
  UINT64 s=0;
  while ((p >> 63) == 0) { p <<= 1; ++s; }
  recint x;
  x.s  = s;
  x.v  = _reciprocal21(p);
  x.d0 = p;
  x.d1 = p;
  return x;
}

inline recint recip2(const UINT64 *p)
{
  UINT64 p0=p[0], p1=p[1];
  if (p1==0) DIV_EXCEPTION();

  UINT64 s=0;
  UINT64 t=p1;
  while ((t >> 63) == 0) { t <<= 1; ++s; }

  if (s) {
    p1 = (p1 << s) | (p0 >> (64-s));
    p0 = (p0 << s);
  }

  recint x;
  x.s  = s;
  x.v  = _reciprocal32(p0,p1);
  x.d0 = p0;
  x.d1 = p1;
  return x;
}

static inline void _rec21(UINT64 u0, UINT64 u1, UINT64 s, UINT64 v, UINT64 d, UINT64 *q, UINT64 *r)
{
  if (s) {
    u1 = (u1 << s) | (u0 >> (64-s));
    u0 = (u0 << s);
  }
  UINT64 q0,q1,r0;
  MUL211(q0,q1, v,u1);
  ADD22(q0,q1,u0,u1+1);
  r0 = u0 - q1*d;
  if (r0 > q0) { q1--; r0 += d; }
  if (r0 >= d) { q1++; r0 -= d; }
  q[0]=q1;
  r[0]=r0 >> s;
}

static inline UINT64 _rec21r(UINT64 u0, UINT64 u1, UINT64 v, UINT64 d)
{
  UINT64 q0,q1,r0;
  MUL211(q0,q1, v,u1);
  ADD22(q0,q1,u0,u1+1);
  r0 = u0 - q1*d;
  if (r0 > q0) r0 += d;
  if (r0 >= d) r0 -= d;
  return r0;
}

static inline void _rec32(UINT64 u0, UINT64 u1, UINT64 u2, UINT64 s, UINT64 v,
                          UINT64 d0, UINT64 d1, UINT64 *q, UINT64 *r)
{
  if (s) {
    u2 = (u2 << s) | (u1 >> (64-s));
    u1 = (u1 << s) | (u0 >> (64-s));
    u0 = (u0 << s);
  }

  UINT64 q0,q1,r0,r1,t0,t1;
  MUL211(q0,q1, v,u2);
  ADD22(q0,q1,u1,u2);

  r0 = u0;
  r1 = u1 - q1*d1;

  MUL211(t0,t1,d0,q1);
  SUB22(r0,r1,d0,d1);
  SUB22(r0,r1,t0,t1);

  q1++;
  if (r1 >= q0) {
    q1--;
    ADD22(r0,r1,d0,d1);
  }
  if (r1 > d1 || (r1==d1 && r0 >= d0)) {
    q1++;
    SUB22(r0,r1,d0,d1);
  }

  if (s) {
    r0 = (r0 >> s) | (r1 << (64-s));
    r1 = (r1 >> s);
  }

  q[0]=q1;
  r[0]=r0; r[1]=r1;
}

static inline void _rec42(UINT64 u0, UINT64 u1, UINT64 u2, UINT64 u3, UINT64 s, UINT64 v,
                          UINT64 d0, UINT64 d1, UINT64 *q, UINT64 *r)
{
  UINT64 q0,q1,R[3];
  _rec32(u1,u2,u3,s,v,d0,d1,&q1,R+1);
  _rec32(u0,R[1],R[2],s,v,d0,d1,&q0,R);
  q[0]=q0; q[1]=q1;
  r[0]=R[0]; r[1]=R[1];
}

inline void rec128(UINT64 *x, recint v)
{
  _rec21(x[0],x[1], v.s, v.v, v.d0, &x[0], &x[1]);
}

inline UINT64 mulrec64(UINT64 a, UINT64 b, recint v)
{
  UINT64 u0,u1,r;
  b <<= v.s;
  MUL211(u0,u1,a,b);
  r = _rec21r(u0,u1,v.v,v.d0);
  return r >> v.s;
}

inline void rec256(UINT64 *x, recint v)
{
  UINT64 q[2], r[2];
  _rec42(x[0],x[1],x[2],x[3], v.s, v.v, v.d0, v.d1, q, r);
  x[0]=q[0]; x[1]=q[1];
  x[2]=r[0]; x[3]=r[1];
}

inline void mulrec128(UINT64 *x, const UINT64 *y, const UINT64 *z, recint v)
{
  UINT64 t[4], q[2], r[2];
  MUL422(t[0],t[1],t[2],t[3],y[0],y[1],z[0],z[1]);
  _rec42(t[0],t[1],t[2],t[3], v.s, v.v, v.d0, v.d1, q, r);
  x[0]=r[0]; x[1]=r[1];
}

/* ------------------- __int128 extensions ------------------- */
#if INT128_EXTENDED

inline void mul256i(UINT64 *x, UINT128 a, UINT128 b)
{
  UINT64 A[2], B[2];
  std::memcpy(A, &a, sizeof(a));
  std::memcpy(B, &b, sizeof(b));
  MUL422(x[0],x[1],x[2],x[3],A[0],A[1],B[0],B[1]);
}

inline void div256i(UINT64 *x, UINT128 a, UINT128 *q, UINT128 *r)
{
  UINT64 A[2], Q[2], R[2];
  std::memcpy(A, &a, sizeof(a));
  _div42(x, A, Q, R);

  UINT128 qq = (UINT128(Q[1]) << 64) | Q[0];
  UINT128 rr = (UINT128(R[1]) << 64) | R[0];
  *q = qq; *r = rr;
}

inline UINT128 mulmod128i(UINT128 a, UINT128 b, UINT128 p)
{
  UINT64 A[2], B[2], P[2], Q[2] = {0,0}, R[2] = {0,0}, x[4] = {0,0,0,0};
  std::memcpy(A, &a, sizeof(a));
  std::memcpy(B, &b, sizeof(b));
  std::memcpy(P, &p, sizeof(p));
  MUL422(x[0],x[1],x[2],x[3],A[0],A[1],B[0],B[1]);
  _div42(x, P, Q, R);
  return (UINT128(R[1]) << 64) | R[0];
}

inline UINT128 mulrec128i(UINT128 a, UINT128 b, recint v)
{
  UINT64 A[2], B[2], Q[2] = {0,0}, R[2] = {0,0}, x[4] = {0,0,0,0};
  std::memcpy(A, &a, sizeof(a));
  std::memcpy(B, &b, sizeof(b));
  MUL422(x[0],x[1],x[2],x[3],A[0],A[1],B[0],B[1]);
  _rec42(x[0],x[1],x[2],x[3], v.s, v.v, v.d0, v.d1, Q, R);
  return (UINT128(R[1]) << 64) | R[0];
}

#endif

/* -------- optional cleanup to avoid macro pollution -------- */
#if !INT128_MACROS
  #undef __ADD
  #undef __ADC
  #undef __SUB
  #undef __SBB
  #undef __MUL
  #undef __DIV
  #undef ADD22
  #undef ADD33
  #undef ADD44
  #undef SUB22
  #undef SUB33
  #undef SUB44
  #undef MUL211
  #undef DIV21H
  #undef FMA211
  #undef FMS211
  #undef MUL321
  #undef FMA321
  #undef FMS321
  #undef MUL422
  #undef FMA422
  #undef FMS422
#endif
