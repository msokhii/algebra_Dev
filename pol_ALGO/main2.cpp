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
// long long cpuMUL=0;
int main(){

    LONG p=9223372036854775783;
    recint P=recip1(p);
    int degN=64;
    int degD=64;
    const int ITER=100000;
    const int STEP=2;

    ofstream logFile("benchMark.txt");
    logFile<<"PRIME -> "<<p<<"\n";
    logFile<<"CALLS -> "<<ITER<<"\n";
    logFile<<left
        <<setw(10)<<"STEP"
        <<setw(10)<<"degN"
        <<setw(10)<<"degD"
        <<setw(18)<<"avgTimeNewton"
        <<setw(18)<<"avgTimeRR"
        <<setw(14)<<"mulsNewton"
        <<setw(14)<<"mulsRR"
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
        int degU=newtonInterp(x.data(),yCopy.data(),m,p,P);
        GLOBALMUL=0;
        auto start=chrono::steady_clock::now();
        for(int i=0;i<ITER;i++){
            copy(y.begin(),y.end(),yCopy.begin());
            int degU=newtonInterp(x.data(),yCopy.data(),m,p,P);
        };
        auto stop=chrono::steady_clock::now();
        auto cpStart=chrono::steady_clock::now();
        for(int i=0;i<ITER;i++){
            copy(y.begin(),y.end(),yCopy.begin());
        }
        auto cpStop=chrono::steady_clock::now();
        double cpTotal=chrono::duration<double,std::micro>(cpStop-cpStart).count();
        double total=chrono::duration<double,std::micro>(stop-start).count();
        double avgTimeCp=cpTotal/ITER;
        double avgTimeNewton=(total/ITER)-avgTimeCp;
        long long newtonMuls=GLOBALMUL/ITER;
        vector<LONG>M(m+1,0);
        M[0]=1;
        int degM=mkM(M,x,p);

        RatReconFastWS W(degM);
        vector<LONG>rOut(m,0);
        vector<LONG>tOut(m,0);
        int degR=-1;
        int degT=-1;
        int flag=-999;
        GLOBALMUL=0;
        auto start2=chrono::steady_clock::now();
        for(int k=0;k<ITER;k++){
            flag=ratReconFastKernelWS(M,y,m,degU,
            degN,degD,p,W,rOut.data(),degR,tOut.data(),degT,P);
        }
        auto stop2=chrono::steady_clock::now();
        long long RRmuls=GLOBALMUL/ITER;
        //long long cpuMULRR=cpuMUL/ITER;
        double total2=chrono::duration<double,std::micro>(stop2-start2).count();
        double avgTimeRR=total2/ITER;
        
        logFile<<left<<
                 setw(10)<<step<<
                 setw(10)<<degN<<
                 setw(10)<<degD<<
                 setw(18)<<avgTimeNewton<<
                 setw(18)<<avgTimeRR<<
                 setw(14)<<newtonMuls<<
                 setw(14)<<RRmuls<<"\n";

        degN*=2;
        degD*=2;
    }
    logFile.close();
    return 0;
}
