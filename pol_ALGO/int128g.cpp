#include "int128g.hpp"

volatile int __int128zero = 0;
volatile int __int128junk = 0;
INT64 CNTR = 0;

/* ------------------- functions ------------------- */

void add128(UINT64 *x, UINT64 *y)
{
    ADD22(x[0], x[1], y[0], y[1]);
}

void sub128(UINT64 *x, UINT64 *y)
{
    SUB22(x[0], x[1], y[0], y[1]);
}

void mul128(UINT64 *x, UINT64 a, UINT64 b)
{
    MUL211(x[0], x[1], a, b);
}

void fma128(UINT64 *x, UINT64 a, UINT64 b)
{
    FMA211(x[0], x[1], a, b);
}

void fms128(UINT64 *x, UINT64 a, UINT64 b)
{
    FMS211(x[0], x[1], a, b);
}

void add256(UINT64 *x, UINT64 *y)
{
    ADD44(x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3]);
}

void sub256(UINT64 *x, UINT64 *y)
{
    SUB44(x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3]);
}

void mul256(UINT64 *x, UINT64 *y, UINT64 *z)
{
    MUL422(x[0], x[1], x[2], x[3], y[0], y[1], z[0], z[1]);
}

void fma256(UINT64 *x, UINT64 *y, UINT64 *z)
{
    FMA422(x[0], x[1], x[2], x[3], y[0], y[1], z[0], z[1]);
}

void fms256(UINT64 *x, UINT64 *y, UINT64 *z)
{
    FMS422(x[0], x[1], x[2], x[3], y[0], y[1], z[0], z[1]);
}

/* ------------------- division ------------------- */

/* 128/64 bit division from Hacker's Delight */
static void _div21(UINT64 u0, UINT64 u1, UINT64 v, UINT64 *q, UINT64 *r)
{
#ifdef DIV21H
    DIV21H(u0, u1, v);
    *q = u0;
    *r = u1;
#else
    UINT64 b, un1, un0, vn1, vn0, q1, q0, un32, un21, un10, rhat;
    INT64 s;

    b = (UINT64)1 << 32;

    /* quotient overflow case */
    if (u1 > v) DIV_EXCEPTION();
    for (s = 0; !(v >> 63); s++) {
        v = v << 1;
    }
    vn1  = v >> 32;
    vn0  = v & 0xFFFFFFFFULL;
    un32 = (u1 << s) | ((u0 >> (64 - s)) & ((UINT64)((-s) >> 63)));
    un10 = u0 << s;
    un1  = un10 >> 32;
    un0  = un10 & 0xFFFFFFFFULL;

    q1   = un32 / vn1;
    rhat = un32 - q1 * vn1;
again1:
    if (q1 >= b || q1 * vn0 > b * rhat + un1) {
        q1 = q1 - 1;
        rhat = rhat + vn1;
        if (rhat < b) goto again1;
    }
    un21 = un32 * b + un1 - q1 * v;
    q0   = un21 / vn1;
    rhat = un21 - q0 * vn1;
again2:
    if (q0 >= b || q0 * vn0 > b * rhat + un0) {
        q0 = q0 - 1;
        rhat = rhat + vn1;
        if (rhat < b) goto again2;
    }
    *r = (un21 * b + un0 - q0 * v) >> s;
    *q = q1 * b + q0;
#endif
}

/* [x0:x1] = [x / a : x % a] */
void div128(UINT64 *x, UINT64 a)
{
    _div21(x[0], x[1], a, &x[0], &x[1]);
}

/* return a*b mod p */
UINT64 mulmod64(UINT64 a, UINT64 b, UINT64 p)
{
    UINT64 x0, x1;
    MUL211(x0, x1, a, b);
    _div21(x0, x1, p, &x0, &x1);
    return x1;
}

/* multiword division from Hacker's Delight */
/* requires that n >= m and v[n-1] non-zero */
static void _divmnu(UINT32 *u, UINT32 *v, UINT32 *q, UINT32 *r, INT32 m, INT32 n)
{
    UINT64 b = 4294967296ULL;
    UINT32 un[m + 1], vn[n];
    UINT64 qhat, rhat, p;
    INT64 t, k;
    INT32 s, i, j;

    while (m && u[m - 1] == 0) m--;
    while (n && v[n - 1] == 0) n--;
    if (n == 0) DIV_EXCEPTION();
    if (n == 1) {
        /* single digit divisor */
        for (k = 0, j = m - 1; j >= 0; j--) {
            q[j] = (k * b + u[j]) / v[0];
            k    = (k * b + u[j]) - (UINT64)q[j] * v[0];
        }
        r[0] = (UINT32)k;
        return;
    }

    /* normalize divisor to set leading bit */
    for (s = 0, t = v[n - 1]; !(t >> 31); t <<= 1, s++) {}
    for (i = n - 1; i > 0; i--) {
        vn[i] = (v[i] << s) | ((UINT64)v[i - 1] >> (32 - s));
    }
    vn[0] = v[0] << s;

    /* shift up dividend */
    un[m] = (UINT64)u[m - 1] >> (32 - s);
    for (i = m - 1; i > 0; i--) {
        un[i] = (u[i] << s) | ((UINT64)u[i - 1] >> (32 - s));
    }
    un[0] = u[0] << s;

    /* multiply and subtract */
    for (j = m - n; j >= 0; j--) {
        /* compute next quotient */
        p    = un[j + n] * b + un[j + n - 1];
        qhat = p / vn[n - 1];
        rhat = p % vn[n - 1];
again:
        if (qhat >= b || qhat * vn[n - 2] > b * rhat + un[j + n - 2]) {
            qhat--;
            rhat += vn[n - 1];
            if (rhat < b) goto again;
        }
        for (k = 0, i = 0; i < n; i++) {
            p      = qhat * vn[i];
            t      = (INT64)un[i + j] - k - (INT64)(p & 0xFFFFFFFFULL);
            un[i + j] = (UINT32)t;
            k      = (p >> 32) - (UINT64)(t >> 32);
        }
        t       = (INT64)un[j + n] - k;
        un[j + n] = (UINT32)t;
        q[j]    = (UINT32)qhat;
        if (t < 0) {
            q[j] = q[j] - 1;
            for (k = 0, i = 0; i < n; i++) {
                t       = (UINT64)un[i + j] + vn[i] + k;
                un[i + j] = (UINT32)t;
                k       = (UINT64)t >> 32;
            }
            un[j + n] += (UINT32)k;
        }
    }

    /* construct remainder */
    for (i = 0; i < n - 1; i++) {
        r[i] = (un[i] >> s) | ((UINT64)un[i + 1] << (32 - s));
    }
    r[n - 1] = un[n - 1] >> s;
}

/* 256/128 bit long division */
/* quotient/remainder of x/y */
static void _div42(UINT64 *x, UINT64 *y, UINT64 *q, UINT64 *r)
{
    UINT32 R[10] = {0};
    UINT32 Q[10] = {0};

    /* quotient overflow */
    if (x[3] > y[1] || (x[3] == 0 && y[1] == 0 && x[2] > y[0])) {
        DIV_EXCEPTION();
    }
    /* divisor greater than dividend */
    if (x[3] == 0 && x[2] == 0 && (x[1] < y[1] || (x[1] == y[1] && x[0] < y[0]))) {
        q[0] = 0;
        q[1] = 0;
        r[0] = x[0];
        r[1] = x[1];
        return;
    }
    _divmnu((UINT32 *)x, (UINT32 *)y, Q, R, 8, 4);
    r[0] = ((UINT64)R[1] << 32) | R[0];
    r[1] = ((UINT64)R[3] << 32) | R[2];
    q[0] = ((UINT64)Q[1] << 32) | Q[0];
    q[1] = ((UINT64)Q[3] << 32) | Q[2];
}

/* [x01:x23] = [x / a : x % a] */
void div256(UINT64 *x, UINT64 *y)
{
    UINT64 q[2], r[2];
    _div42(x, y, q, r);
    x[0] = q[0];
    x[1] = q[1];
    x[2] = r[0];
    x[3] = r[1];
}

/* [x0:x1] = [y0:y1] * [z0:z1] mod [p0:p1] */
void mulmod128(UINT64 *x, UINT64 *y, UINT64 *z, UINT64 *p)
{
    UINT64 t[4], q[2], r[2];
    MUL422(t[0], t[1], t[2], t[3], y[0], y[1], z[0], z[1]);
    _div42(t, p, q, r);
    x[0] = r[0];
    x[1] = r[1];
}

/* --------------- inverse division --------------- */

/* for 2/1 division */
static UINT64 _reciprocal21(UINT64 p)
{
    UINT64 u0, u1, q, r;
    u1 = -(p + 1);
    u0 = (UINT64)-1;
    _div21(u0, u1, p, &q, &r);
    return q;
}

/* for 3/2 division */
static UINT64 _reciprocal32(UINT64 d0, UINT64 d1)
{
    UINT64 t0, t1, p, v;
    v = _reciprocal21(d1);
    p = d1 * v + d0;
    if (p < d0) {
        v--;
        if (p >= d1) {
            v--;
            p -= d1;
        }
        p -= d1;
    }
    MUL211(t0, t1, v, d0);
    p += t1;
    if (p < t1) {
        v--;
        if (p >= d1 || (p == d1 && t0 >= d0)) {
            v--;
        }
    }
    return v;
}

/* 64-bit reciprocal */
recint recip1(UINT64 p)
{
    UINT64 s;
    recint x;
    if (p == 0) DIV_EXCEPTION();
    for (s = 0; !(p >> 63); p <<= 1, s++) {}
    x.s  = s;
    x.v  = _reciprocal21(p);
    x.d0 = p;
    x.d1 = p;
    return x;
}

/* 128-bit reciprocal */
recint recip2(UINT64 *p)
{
    UINT64 p0, p1, t, s;
    recint x;
    p0 = p[0];
    p1 = p[1];
    if (p1 == 0) DIV_EXCEPTION();
    for (s = 0, t = p1; !(t >> 63); t <<= 1, s++) {}
    if (s) {
        p1 = (p1 << s) | (p0 >> (64 - s));
        p0 = (p0 << s);
    }
    x.s  = s;
    x.v  = _reciprocal32(p0, p1);
    x.d0 = p0;
    x.d1 = p1;
    return x;
}

/* 128/64 bit division */
static void _rec21(UINT64 u0, UINT64 u1, UINT64 s, UINT64 v, UINT64 d, UINT64 *q, UINT64 *r)
{
    UINT64 q0, q1, r0;
    if (s) {
        u1 = (u1 << s) | (u0 >> (64 - s));
        u0 = (u0 << s);
    }
    MUL211(q0, q1, v, u1);
    ADD22(q0, q1, u0, u1 + 1);
    r0 = u0 - q1 * d;
    if (r0 > q0) {
        q1--;
        r0 += d;
    }
    if (r0 >= d) {
        q1++;
        r0 -= d;
    }
    q[0] = q1;
    r[0] = r0 >> s;
}

/* 128/64 bit remainder, [u0:u1] already shifted up */
static UINT64 _rec21r(UINT64 u0, UINT64 u1, UINT64 v, UINT64 d)
{
    UINT64 q0, q1, r0;
    MUL211(q0, q1, v, u1);
    ADD22(q0, q1, u0, u1 + 1);
    r0 = u0 - q1 * d;
    if (r0 > q0) {
        r0 += d;
    }
    if (r0 >= d) {
        r0 -= d;
    }
    return r0;
}

/* 192/128 bit division */
static void _rec32(UINT64 u0, UINT64 u1, UINT64 u2, UINT64 s, UINT64 v,
                   UINT64 d0, UINT64 d1, UINT64 *q, UINT64 *r)
{
    UINT64 q0, q1, r0, r1, t0, t1;
    if (s) {
        u2 = (u2 << s) | (u1 >> (64 - s));
        u1 = (u1 << s) | (u0 >> (64 - s));
        u0 = (u0 << s);
    }
    MUL211(q0, q1, v, u2);
    ADD22(q0, q1, u1, u2);
    r0 = u0;
    r1 = u1 - q1 * d1;
    MUL211(t0, t1, d0, q1);
    SUB22(r0, r1, d0, d1);
    SUB22(r0, r1, t0, t1);
    q1++;
    if (r1 >= q0) {
        q1--;
        ADD22(r0, r1, d0, d1);
    }
    if (r1 > d1 || ((r1 == d1) && (r0 >= d0))) {
        q1++;
        SUB22(r0, r1, d0, d1);
    }
    if (s) {
        r0 = (r0 >> s) | (r1 << (64 - s));
        r1 = (r1 >> s);
    }
    q[0] = q1;
    r[0] = r0;
    r[1] = r1;
}

/* 256/128 bit division */
static void _rec42(UINT64 u0, UINT64 u1, UINT64 u2, UINT64 u3, UINT64 s, UINT64 v,
                   UINT64 d0, UINT64 d1, UINT64 *q, UINT64 *r)
{
    UINT64 q0, q1, R[3];
    _rec32(u1, u2, u3, s, v, d0, d1, &q1, R + 1);
    _rec32(u0, R[1], R[2], s, v, d0, d1, &q0, R);
    q[0] = q0;
    q[1] = q1;
    r[0] = R[0];
    r[1] = R[1];
}

/* 128/64 bit division */
void rec128(UINT64 *x, recint v)
{
    _rec21(x[0], x[1], v.s, v.v, v.d0, &x[0], &x[1]);
}

/* 64 x 64 bit mulmod */
UINT64 mulrec64(UINT64 a, UINT64 b, recint v)
{
    UINT64 u0, u1, r;
    b = b << v.s;
    MUL211(u0, u1, a, b);
    r = _rec21r(u0, u1, v.v, v.d0);
    return r >> v.s;
}

/* 256/128 bit division */
void rec256(UINT64 *x, recint v)
{
    UINT64 q[2], r[2];
    _rec42(x[0], x[1], x[2], x[3], v.s, v.v, v.d0, v.d1, q, r);
    x[0] = q[0];
    x[1] = q[1];
    x[2] = r[0];
    x[3] = r[1];
}

/* 128 x 128 bit mulmod */
void mulrec128(UINT64 *x, UINT64 *y, UINT64 *z, recint v)
{
    UINT64 t[4], q[2], r[2];
    MUL422(t[0], t[1], t[2], t[3], y[0], y[1], z[0], z[1]);
    _rec42(t[0], t[1], t[2], t[3], v.s, v.v, v.d0, v.d1, q, r);
    x[0] = r[0];
    x[1] = r[1];
}

/* --------------- 128-bit integers --------------- */
#if INT128_EXTENDED

void mul256i(UINT64 *x, UINT128 a, UINT128 b)
{
    UINT64 A[2], B[2];
    A[0] = *((UINT64 *)(&a) + 0);
    A[1] = *((UINT64 *)(&a) + 1);
    B[0] = *((UINT64 *)(&b) + 0);
    B[1] = *((UINT64 *)(&b) + 1);
    MUL422(x[0], x[1], x[2], x[3], A[0], A[1], B[0], B[1]);
}

void div256i(UINT64 *x, UINT128 a, UINT128 *q, UINT128 *r)
{
    UINT128 q1, r1;
    UINT64 A[2], Q[2], R[2];
    A[0] = *((UINT64 *)(&a) + 0);
    A[1] = *((UINT64 *)(&a) + 1);
    _div42(x, A, Q, R);
    q1 = Q[1];
    q1 = (q1 << 64) | Q[0];
    r1 = R[1];
    r1 = (r1 << 64) | R[0];
    *q = q1;
    *r = r1;
}

UINT128 mulmod128i(UINT128 a, UINT128 b, UINT128 p)
{
    UINT128 r;
    UINT64 A[2], B[2], P[2], Q[2] = {0, 0}, R[2] = {0, 0}, x[4] = {0, 0, 0, 0};
    A[0] = *((UINT64 *)(&a) + 0);
    A[1] = *((UINT64 *)(&a) + 1);
    B[0] = *((UINT64 *)(&b) + 0);
    B[1] = *((UINT64 *)(&b) + 1);
    P[0] = *((UINT64 *)(&p) + 0);
    P[1] = *((UINT64 *)(&p) + 1);
    MUL422(x[0], x[1], x[2], x[3], A[0], A[1], B[0], B[1]);
    _div42(x, P, Q, R);
    r = R[1];
    r = (r << 64) | R[0];
    return r;
}

UINT128 mulrec128i(UINT128 a, UINT128 b, recint v)
{
    UINT128 r;
    UINT64 A[2], B[2], Q[2] = {0, 0}, R[2] = {0, 0}, x[4] = {0, 0, 0, 0};
    A[0] = *((UINT64 *)(&a) + 0);
    A[1] = *((UINT64 *)(&a) + 1);
    B[0] = *((UINT64 *)(&b) + 0);
    B[1] = *((UINT64 *)(&b) + 1);
    MUL422(x[0], x[1], x[2], x[3], A[0], A[1], B[0], B[1]);
    _rec42(x[0], x[1], x[2], x[3], v.s, v.v, v.d0, v.d1, Q, R);
    r = R[1];
    r = (r << 64) | R[0];
    return r;
}

#endif