#include"polyMath.h"
#include"integerMath.h"
#include"interpAlgo.h"
#include"rfrWrapper.h"
#include<vector>
#include<iostream>
#include<cstdint>
#include<time.h>
#include<chrono>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

using namespace std; 

int mkM(vector<LONG>&m,const vector<LONG> &xs,const LONG p){    
    int degM=0;
    std::vector<LONG>linF(2,0);
    linF[1]=1;
    for(int i=0;i<xs.size();i++){
        linF[0]=(xs[i]==0?0:p-xs[i]);
        degM=pMULIP64(m.data(),linF.data(),degM,1,p);
    };
    return degM;
}

long long GLOBALMUL=0;
long long GLOBALMUL64=0;
long long GLOBALCPUMUL=0;
int main(){

    LONG p=9223372036854775783;
    recint P=recip1(p);
    int degN=5;
    int degD=5;
    const int ITER=100000;
    const int STEP=5;

    ofstream logFile("benchMark.txt");
    logFile<<"PRIME -> "<<p<<"\n";
    logFile<<"CALLS -> "<<ITER<<"\n";
    logFile<<left
        <<setw(10)<<"STEP"
        <<setw(10)<<"degN"
        <<setw(10)<<"degD"
        <<setw(28)<<"avgTimeNewton(mulRec)"
<<<<<<< HEAD
        <<setw(28)<<"avgTimeNewton(mulSub)"
        <<setw(28)<<"avgTimeRR(No CPU ~ mulRec)"
        <<setw(28)<<"avgTimeRR(CPU ~ mulSub)"
=======
        <<setw(28)<<"avgTimeNewton(mul64)"
        <<setw(28)<<"avgTimeRR(No CPU+mulRec)"
        <<setw(28)<<"avgTimeRR(CPU+mul64)"
>>>>>>> 391ff8f547858484a25f15ecfd751e1cd5abd64f
        <<setw(22)<<"mulsNewton(mulRec)"
        <<setw(22)<<"mulsNewton(mulSub)"
        <<setw(18)<<"mulsRR"
        <<setw(18)<<"mulsRRNorm"
        << "\n";

    for(int step=1;step<STEP;step++){
        vector<LONG>n(degN+1,0);
        vector<LONG>d(degD+1,0);

        for(int i=0;i<degN+1;i++){
            LONG temp=rand64s(p);
            while(temp==0){
                temp=rand64s(p);
            }
            n[i]=temp;
        }
        for(int j=0;j<degD+1;j++){
            LONG temp=rand64s(p);
            while(temp==0){
                temp=rand64s(p);
            }
            d[j]=temp;
        }

        if(d[degD]!=1){
            LONG invTerm;
            invTerm=modinv64b(d[degD],p);
            for(int i=0;i<=degD;i++){
                d[i]=mul64b(invTerm,d[i],p);
            }
            for(int j=0;j<=degN;j++){
                n[j]=mul64b(invTerm,n[j],p);
            }
        }

        vector<LONG>nCopy=n;
        vector<LONG>dCopy=d;

        int m=degN+degD+1;
        vector<LONG>x(m,0);
        for(int i=0;i<m;i++){
            x[i]=i+1;
        }

        vector<LONG>y(m,0);
        for(int i=0;i<m;i++){
            LONG denEval=pEVAL64(d.data(),degD,x[i],p);
            if(denEval==0){
                return -1;
            }
            LONG numEval=pEVAL64(n.data(),degN,x[i],p);
            y[i]=mul64b(numEval,modinv64b(denEval,p),p);
        }

        vector<LONG>yCopy(m,0);
        copy(y.begin(),y.end(),yCopy.begin());
        int degU=newtonInterpMulRec(x.data(),yCopy.data(),m,p,P);
        
        // Timer for newton interpolation using mulrec64 routine.
        GLOBALMUL=0;
        auto start=chrono::steady_clock::now();
        for(int i=0;i<ITER;i++){
            copy(y.begin(),y.end(),yCopy.begin());
            int degU=newtonInterpMulRec(x.data(),yCopy.data(),m,p,P);
        };
        auto stop=chrono::steady_clock::now();
        
        // Timer for newton interpolation using mul64b routine.
        GLOBALMUL64=0;
        auto newton2Start=chrono::steady_clock::now();
        for(int i=0;i<ITER;i++){
            copy(y.begin(),y.end(),yCopy.begin());
            int degU=newtonInterpMulNormal(x.data(),yCopy.data(),m,p);
        }
        auto newton2Stop=chrono::steady_clock::now();
        
        // Timer for copying y into y0 for newton interpolation.
        auto cpStart=chrono::steady_clock::now();
        for(int i=0;i<ITER;i++){
            copy(y.begin(),y.end(),yCopy.begin());
        }
        auto cpStop=chrono::steady_clock::now();
        double cpTotal=chrono::duration<double,std::micro>(cpStop-cpStart).count();
        double total=chrono::duration<double,std::micro>(stop-start).count();
        double newton2Total=chrono::duration<double,std::micro>(newton2Stop-newton2Start).count();
        double avgTimeCp=cpTotal/ITER;
        double avgTimeNewton=(total/ITER)-avgTimeCp;
        double avgTimeNewton2=(newton2Total/ITER)-avgTimeCp;
        long long newtonMuls=GLOBALMUL/ITER;
        long long newtonMuls2=GLOBALMUL64/ITER;
        vector<LONG>M(m+1,0);
        M[0]=1;
        int degM=mkM(M,x,p);

        RatReconFastWS W(degM);
        RatReconFastWS W2(degM);
        vector<LONG>rOut(m,0);
        vector<LONG>tOut(m,0);
        vector<LONG>rOut2(m,0);
        vector<LONG>tOut2(m,0);
        int degR=-1;
        int degT=-1;
        int flag=-999;
        int degR2=-1;
        int degT2=-1;
        int flag2=-999;
        int degUCP=degU;
        int degNCP=degN;
        int degDCP=degD;
        int mCP=m;
        vector<LONG> MCP=M;
        vector<LONG> yCP2=y;
        GLOBALMUL=0;
        auto start2=chrono::steady_clock::now();
        for(int k=0;k<ITER;k++){
            flag=ratReconFastKernelWS(M,y,m,degU,
            degN,degD,p,W,rOut.data(),degR,tOut.data(),degT,P);
        }
        auto stop2=chrono::steady_clock::now();
<<<<<<< HEAD
        GLOBALMUL64=0;
        GLOBALCPUMUL=0;
=======
        GLOBALCPUMUL=0;
        GLOBALMUL64=0;
>>>>>>> 391ff8f547858484a25f15ecfd751e1cd5abd64f
        auto rrNormStart=chrono::steady_clock::now();
        for(int k=0;k<ITER;k++){
            flag=ratReconNormal(MCP,yCP2,mCP,degUCP,
            degNCP,degDCP,p,W2,rOut2.data(),degR2,tOut2.data(),degT2);
        }
        auto rrNormStop=chrono::steady_clock::now();

        long long RRmuls=GLOBALMUL/ITER;
        long long cpuMULRR=GLOBALCPUMUL/ITER;
        long long totalRRMuls=RRmuls;
        long long rrMulsNorm=(GLOBALMUL64/ITER)+cpuMULRR;
        double total2=chrono::duration<double,std::micro>(stop2-start2).count();
        double rrNormTotal=chrono::duration<double,std::micro>(rrNormStop-rrNormStart).count();
        double avgTimeRR=total2/ITER;
        double avgTimeRRNorm=rrNormTotal/ITER;

        logFile<<left<<
                 setw(10)<<step<<
                 setw(10)<<degN<<
                 setw(10)<<degD<<
                 setw(28)<<avgTimeNewton<<
                 setw(28)<<avgTimeNewton2<<
                 setw(28)<<avgTimeRR<<
                 setw(28)<<avgTimeRRNorm<<
                 setw(22)<<newtonMuls<<
                 setw(22)<<newtonMuls2<<
                 setw(18)<<totalRRMuls<<
                 setw(18)<<rrMulsNorm<<
                 "\n";
                  
        degN*=2;
        degD*=2;
    }
    logFile.close();
    return 0;
}
