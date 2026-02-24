/* 

VERSION 1: 

#include <iostream>
#include <vector>
#include <chrono>
#include <stdexcept>
#include "polyMath.h"
#include "integerMath.h"

using namespace std;

void vSOLVER64(
    const vector<LONG>& m,
    const vector<LONG>& y,
    int n,
    vector<LONG>& a,
    vector<LONG>& M,
    int shift,
    LONG p
) {
    if ((int)m.size() < n || (int)y.size() < n) {
        throw std::invalid_argument("m and y must have at least n entries");
    }
    if (n <= 0) {
        a.clear();
        M.clear();
        return;
    }

    auto normp = [p](LONG x) -> LONG {
        x %= p;
        if (x < 0) x += p;
        return x;
    };

    // Normalize inputs
    vector<LONG> mm(n), yy(n);
    for (int i = 0; i < n; ++i) {
        mm[i] = normp(m[i]);
        yy[i] = normp(y[i]);
    }

    // Check duplicate nodes mod p
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (mm[i] == mm[j]) {
                cout << "roots are not distinct (duplicate nodes mod p)\n";
                a.assign(n, 0);
                M.assign(n + 1, 0);
                return;
            }
        }
    }

    a.assign(n, 0);
    M.assign(n + 1, 0);

    vector<LONG> Q(n, 0);
    vector<LONG> Prod(n + 1, 0);
    vector<LONG> T(n + 1, 0);
    vector<LONG> A(2, 0);

    // Safe in-place multiply by (x - r) for low->high coeffs
    auto mulByLinearSafe = [&](vector<LONG>& P, int degP, LONG r) {
        vector<LONG> R(degP + 2, 0);
        LONG c0 = neg64s(r, p); // constant term = -r mod p

        // (P0 + P1 x + ... + P_deg x^deg) * (c0 + x)
        for (int i = 0; i <= degP; ++i) {
            R[i]   = add64b(R[i],   mul64bASM(P[i], c0, p), p); // * c0
            R[i+1] = add64b(R[i+1], P[i], p);                   // * x
        }

        for (int i = 0; i <= degP + 1; ++i) P[i] = R[i];
    };

    Prod[0] = 1;

    auto t0 = chrono::high_resolution_clock::now();

    // Build Prod(x) = Π (x - mm[i]) safely
    for (int i = 0; i < n; ++i) {
        mulByLinearSafe(Prod, i, mm[i]);
    }

    auto t1 = chrono::high_resolution_clock::now();

    // Return product polynomial in M
    M = Prod;

    cout << "Mul time = "
         << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count()
         << " ms\n";

    // Sanity check: each mm[j] should be a root of Prod
    for (int j = 0; j < n; ++j) {
        LONG chk = evalHORN64(Prod, mm[j], p);
        if (chk != 0) {
            cout << "Product root check failed at j=" << j
                 << "  Prod(m[j])=" << chk << "\n";
            return;
        }
    }

    // Solve
    for (int j = 0; j < n; ++j) {
        A[0] = neg64s(mm[j], p);
        A[1] = 1;

        T = Prod;
        polDIVIP64(T, A, n, 1, p);

        // For your polDIVIP64:
        //   T[0] = remainder, T[1..n] = quotient
        if (T[0] != 0) {
            cout << "Warning: nonzero remainder at j=" << j
                 << " remainder=" << T[0] << "\n";
            return;
        }

        for (int i = 0; i < n; ++i) {
            Q[i] = T[i + 1];
        }

        LONG u = evalHORN64(Q, mm[j], p);
        if (u == 0) {
            cout << "roots are not distinct\n";
            return;
        }

        u = modinv64b(u, p);

        LONG s = 0;
        for (int i = 0; i < n; ++i) {
            s = add64b(s, mul64bASM(Q[i], yy[i], p), p);
        }

        s = mul64bASM(u, s, p);

        if (shift != 0) {
            LONG invmj = modinv64b(mm[j], p);
            LONG pw = powmod64s(invmj, shift, p);
            s = mul64bASM(pw, s, p);
        }

        a[j] = s;
    }

    auto t2 = chrono::high_resolution_clock::now();
    cout << "Solve time = "
         << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count()
         << " ms\n";
}

*/

#include <iostream>
#include <vector>
#include <chrono>
#include <stdexcept>
#include "polyMath.h"
#include "integerMath.h"

using namespace std;

void vSOLVER64(
    const vector<LONG>& m,
    const vector<LONG>& y,
    int n,
    vector<LONG>& a,
    vector<LONG>& M,
    int shift,
    LONG p
) {
    if ((int)m.size() < n || (int)y.size() < n) {
        throw std::invalid_argument("m and y must have at least n entries");
    }
    if (n <= 0) {
        a.clear();
        M.clear();
        return;
    }

    auto normp = [p](LONG x) -> LONG {
        x %= p;
        if (x < 0) x += p;
        return x;
    };

    // Normalize inputs
    vector<LONG> mm(n), yy(n);
    for (int i = 0; i < n; ++i) {
        mm[i] = normp(m[i]);
        yy[i] = normp(y[i]);
    }

    // Check duplicate nodes mod p
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (mm[i] == mm[j]) {
                cout << "roots are not distinct (duplicate nodes mod p)\n";
                a.assign(n, 0);
                M.assign(n + 1, 0);
                return;
            }
        }
    }

    a.assign(n, 0);
    M.assign(n + 1, 0);

    vector<LONG> Q(n, 0);
    vector<LONG> T(n + 1, 0);
    vector<LONG> A(2, 0);

    // Product polynomial Prod(x) = \prod_i (x - mm[i])
    // Keep it dynamic; pMULIP64 will grow it as needed.
    vector<LONG> Prod(1, 1);   // degree 0 polynomial "1"

    auto t0 = chrono::high_resolution_clock::now();

    // Build Prod(x) = Π (x - mm[i]) using your new pMULIP64 (b <- a*b)
    A[1] = 1;
    for (int i = 0; i < n; ++i) {
        A[0] = neg64s(mm[i], p);      // A(x) = x - mm[i]
        pMULIP64VANDER(A, Prod, 1, i, p);   // Prod <- A * Prod
    }

    auto t1 = chrono::high_resolution_clock::now();

    // Return product polynomial in M (size n+1 expected by caller)
    M.assign(n + 1, 0);
    for (int i = 0; i < (int)Prod.size() && i <= n; ++i) {
        M[i] = Prod[i];
    }

    cout << "Mul time = "
         << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count()
         << " ms\n";

    // Sanity check: each mm[j] should be a root of Prod
    for (int j = 0; j < n; ++j) {
        LONG chk = evalHORN64(Prod, mm[j], p);
        if (chk != 0) {
            cout << "Product root check failed at j=" << j
                 << "  Prod(m[j])=" << chk << "\n";
            return;
        }
    }

    // Solve
    for (int j = 0; j < n; ++j) {
        A[0] = neg64s(mm[j], p);
        A[1] = 1;

        T = Prod;                    // copy full product
        polDIVIP64(T, A, n, 1, p);   // divide by (x - mm[j])

        // For your polDIVIP64:
        //   T[0] = remainder, T[1..n] = quotient
        if (T[0] != 0) {
            cout << "Warning: nonzero remainder at j=" << j
                 << " remainder=" << T[0] << "\n";
            return;
        }

        for (int i = 0; i < n; ++i) {
            Q[i] = T[i + 1];
        }

        LONG u = evalHORN64(Q, mm[j], p);
        if (u == 0) {
            cout << "roots are not distinct\n";
            return;
        }

        u = modinv64b(u, p);

        LONG s = 0;
        for (int i = 0; i < n; ++i) {
            s = add64b(s, mul64bASM(Q[i], yy[i], p), p);
        }

        s = mul64bASM(u, s, p);

        if (shift != 0) {
            LONG invmj = modinv64b(mm[j], p);
            LONG pw = powmod64s(invmj, shift, p);
            s = mul64bASM(pw, s, p);
        }

        a[j] = s;
    }

    auto t2 = chrono::high_resolution_clock::now();
    cout << "Solve time = "
         << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count()
         << " ms\n";
}