#include"polyMath.h"
#include"integerMath.h"
#include<vector>
#include<cstdint>
#include"fastVSolve.h"

using namespace std;
using LONG=int64_t;

extern "C" void vSOLVER_C(const LONG* m_in,
                          const LONG* y_in,
                          int n,
                          LONG* a_out,
                          LONG* M_out,
                          int shift,
                          LONG p){
    if (!m_in||!y_in||!a_out||!M_out||n<=0){return;}
    vector<LONG> m(m_in,m_in+n);
    vector<LONG> y(y_in,y_in+n);
    vector<LONG> a(n,0);
    vector<LONG> M(n+1,0); 
    vSOLVER64(m, y, n, a, M, shift, p);
    for (int i=0;i<n;i++){
        a_out[i]=a[i];
    }
    for (int i=0;i<=n;i++){
        M_out[i]=M[i];
    }
}