#include"integerMath.h"
#include"helperF.h"

using namespace std; 

using LONG=int_fast64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

void makeMonic(vector<LONG> &v,int deg,LONG p){
	if(deg<0) return;
	LONG LC=v[deg];
	if(LC==1) return;
	LONG invLC=modinv64b(LC,p);
	for(int i=0;i<=deg;i++){
		v[i]=mul64b(v[i],invLC,p);
	}
}