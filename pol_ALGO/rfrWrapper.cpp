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
    // -----------------------------
    // basic input checks
    // -----------------------------
    if (!M || !U || !nOut || !dOut || !degNOUT || !degDOUT) {
        return -1;
    }
    if (degM < 0 || degU < 0) {
        return -1;
    }
    if (mLen <= 0 || uLen <= 0 || nOutLen <= 0 || dOutLen <= 0) {
        return -1;
    }
    if (degM >= mLen || degU >= uLen) {
        return -1;
    }

    *degNOUT = -1;
    *degDOUT = -1;

    // zero user outputs first
    for (int i = 0; i < nOutLen; ++i) nOut[i] = 0;
    for (int i = 0; i < dOutLen; ++i) dOut[i] = 0;

    // -----------------------------
    // copy input arrays into vectors
    // -----------------------------
    
    vector<LONG> m(M, M + (degM + 1));
    vector<LONG> u(U, U + (degU + 1));

    // -----------------------------
    // workspace size
    // Use something safely large enough for all internal polys.
    // -----------------------------
    int wsSize = std::max(degM, degU) + 1;
    RatReconFastWS W(wsSize);

    // reciprocal structure for mulrec64
    recint P = recip1(p);

    // -----------------------------
    // temporary output buffers
    // do NOT write directly into user buffers yet
    // -----------------------------
    vector<LONG> rTmp(wsSize, 0);
    vector<LONG> tTmp(wsSize, 0);

    int degROut = -1;
    int degTOut = -1;

    // -----------------------------
    // call kernel
    // -----------------------------
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
        printf("HI");
        return rc;
    }

    // -----------------------------
    // sanity checks on kernel output
    // -----------------------------
    if (degROut < 0 || degTOut < 0) {
        *degNOUT = -1;
        *degDOUT = -1;
        return -2;
    }

    if (degROut >= wsSize || degTOut >= wsSize) {
        *degNOUT = -1;
        *degDOUT = -1;
        return -3;
    }

    if (degROut + 1 > nOutLen) {
        *degNOUT = -1;
        *degDOUT = -1;
        return -4;
    }

    if (degTOut + 1 > dOutLen) {
        *degNOUT = -1;
        *degDOUT = -1;
        return -5;
    }

    // -----------------------------
    // copy temp outputs to user buffers
    // -----------------------------
    std::copy_n(rTmp.data(), degROut + 1, nOut);
    std::copy_n(tTmp.data(), degTOut + 1, dOut);

    *degNOUT = degROut;
    *degDOUT = degTOut;

    return 0;
}
