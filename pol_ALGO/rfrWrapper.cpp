#include"polyMath.h"
#include"integerMath.h"
#include<algorithm>
#include<vector>
#include<cstdint> 

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
    // basic input checks
    if (!M || !U || !nOut || !dOut || !degNOUT || !degDOUT) return -1;
    if (degM < 0 || degU < 0) return -1;
    if (degM >= mLen || degU >= uLen) return -1;
    if (nOutLen <= 0 || dOutLen <= 0) return -1;

    // make local vectors from input arrays
    vector<LONG> m(M, M + (degM + 1));
    vector<LONG> u(U, U + (degU + 1));

    // workspace
    RatReconFastWS W(degM+1);

    // output degrees
    int degROut = -1;
    int degTOut = -1;

    // optional: zero outputs first
    for (int i = 0; i < nOutLen; ++i) nOut[i] = 0;
    for (int i = 0; i < dOutLen; ++i) dOut[i] = 0;

    // call kernel
    int rc = ratReconFastKernelWS(m,
                                  u,
                                  degM,
                                  degU,
                                  N,
                                  D,
                                  p,
                                  W,
                                  nOut,
                                  degROut,
                                  dOut,
                                  degTOut);

    if (rc != 0) {
        *degNOUT = -1;
        *degDOUT = -1;
        return rc;
    }

    // bounds check on produced output
    if (degROut + 1 > nOutLen) {
        *degNOUT = -1;
        *degDOUT = -1;
        return -2;
    }
    if (degTOut + 1 > dOutLen) {
        *degNOUT = -1;
        *degDOUT = -1;
        return -3;
    }

    *degNOUT = degROut;
    *degDOUT = degTOut;
    return 0;
}
