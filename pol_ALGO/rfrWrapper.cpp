#include "polyMath.h"
#include "integerMath.h"
#include "int128g.hpp"
#include <algorithm>
#include <vector>
#include <cstdint>
#include <time.h>
#include <chrono>

using namespace std;
using LONG = int64_t;

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
    using clock = std::chrono::steady_clock;

    static long long calls = 0;
    static long long time_m_ns = 0;
    static long long time_u_ns = 0;
    static long long time_r_ns = 0;
    static long long time_t_ns = 0;
    static long long time_struct = 0;
    static long long cp1 = 0;
    static long long cp2 = 0;

    if (!M || !U || !nOut || !dOut || !degNOUT || !degDOUT) return -1;
    if (degM < 0 || degU < 0) return -1;
    if (mLen <= 0 || uLen <= 0 || nOutLen <= 0 || dOutLen <= 0) return -1;
    if (degM >= mLen || degU >= uLen) return -1;

    *degNOUT = -1;
    *degDOUT = -1;

    for (int i = 0; i < nOutLen; ++i) nOut[i] = 0;
    for (int i = 0; i < dOutLen; ++i) dOut[i] = 0;

    int wsSize = std::max(degM, degU) + 1;

    auto t0 = clock::now();
    std::vector<LONG> m(M, M + (degM + 1));
    auto t1 = clock::now();
    time_m_ns += chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    t0 = clock::now();
    std::vector<LONG> u(U, U + (degU + 1));
    t1 = clock::now();
    time_u_ns += chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    t0 = clock::now();
    std::vector<LONG> rTmp(wsSize, 0);
    t1 = clock::now();
    time_r_ns += chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    t0 = clock::now();
    std::vector<LONG> tTmp(wsSize, 0);
    t1 = clock::now();
    time_t_ns += chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    int degROut = -1;
    int degTOut = -1;
    t0 = clock::now();
    RatReconFastWS W(wsSize);
    t1 = clock::now();
    time_struct += chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
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

    calls++;

    if (calls % 1000 == 0) {
        printf("Average vector times: \n");
        printf("m    : %.5f ns\n", (double)time_m_ns / calls);
        printf("u    : %.5f ns\n", (double)time_u_ns / calls);
        printf("rTmp : %.5f ns\n", (double)time_r_ns / calls);
        printf("tTmp : %.5f ns\n", (double)time_t_ns / calls);
        printf("tStruct : %.5f ns\n", (double)time_struct / calls);
        printf("CP1 : %.5f ns\n",(double)cp1/calls);
        printf("CP2: %.5f ns\n",(double)cp2/calls);
        fflush(stdout);
    }

    if (rc != 0) {
        *degNOUT = -1;
        *degDOUT = -1;
        return rc;
    }

    if (degROut < 0 || degTOut < 0) return -2;
    if (degROut >= wsSize || degTOut >= wsSize) return -3;
    if (degROut + 1 > nOutLen) return -4;
    if (degTOut + 1 > dOutLen) return -5;
    t0 = clock::now();
    std::copy_n(rTmp.data(), degROut + 1, nOut);
    t1 = clock::now();
    cp1 += std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
    t0 = clock::now();
    std::copy_n(tTmp.data(), degTOut + 1, dOut);
    t1 = clock::now();
    cp2 += std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();

    *degNOUT = degROut;
    *degDOUT = degTOut;

    return 0;
}
