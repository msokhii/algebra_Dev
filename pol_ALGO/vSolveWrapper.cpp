#include <vector>
#include"polyMath.h"
#include"integerMath.h"
#include"fastVSolve.h"

using std::vector;

// C-compatible exported wrapper (for Maple define_external / C calls)
extern "C" void vSOLVER64C(
    LONG* m_in,   // length n
    LONG* y_in,   // length n
    int   n,
    LONG* a_out,  // length n
    LONG* M_out,  // length n+1
    int   shift,
    LONG  p
) {
    if (!m_in || !y_in || !a_out || !M_out || n <= 0) return;

    // raw arrays -> vectors
    vector<LONG> m(m_in, m_in + n);
    vector<LONG> y(y_in, y_in + n);
    auto normp = [p](LONG x) -> LONG {
    x %= p;
    if (x < 0) x += p;
    return x;
};

for (int i = 0; i < n; ++i) {
    m[i] = normp(m[i]);
    y[i] = normp(y[i]);
}
    vector<LONG> a;
    vector<LONG> M;

    // call your C++ solver
    vSOLVER64(y,m, n, a, M, shift, p);

    // copy results back out
    for (int i = 0; i < n; ++i) {
        a_out[i] = a[i];
    }
    for (int i = 0; i < n + 1; ++i) {
        M_out[i] = M[i];
    }
}