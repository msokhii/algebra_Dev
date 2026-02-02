/* 
Step 1: Raw cpp file. Compile it in cpp to make sure it works. 
Maple expects a C style ABI i.e. can only bind to raw pointers and 
explicit values and return array should be pre-allocated. 

Step 2: Remove the main() and just compile using: 
g++ -O3 -shared -fPIC -o fileName.so fileName.cpp

To check dependencies: ldd fileName.so 

Step 3: Compile using maple command: 
maple fileName.mpl
*/

#include<cstdint>

extern "C" void testAdd64s(const int64_t *A,const int64_t *B,int32_t n,int64_t *C){
    for(int32_t i=0;i<n;i++){
        C[i]=A[i]+B[i];
    }
}