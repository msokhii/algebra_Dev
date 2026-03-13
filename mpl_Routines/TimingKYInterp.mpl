(* First need to compile run.sh and runWrap2.sh and then call maple _FILE_NAME_ *)

restart:

with(NumberTheory):
with(LinearAlgebra):

libV := "/cecm/home/mss59/Desktop/update/pol_ALGO/fastVSolve.so":
libR := "/cecm/home/mss59/Desktop/update/pol_ALGO/rfr.so":

(* Using 0-based indexing *)
mVSOLVE := define_external(
                           'vSOLVER64C',
                           RR::ARRAY(0..nn-1,datatype=integer[8]),
                           LL::ARRAY(0..nn-1,datatype=integer[8]),
                           nn::integer[4],
                           aa::ARRAY(0..nn-1,datatype=integer[8]),
                           XX::ARRAY(0..nn,datatype=integer[8]),
                           shift::integer[4],
                           pp::integer[8],
                           LIB=libV
):

(* Using 0-based indexing *)
mRATRECON := define_external(
                            'ratRECON_C',
                            mLen::integer[4],
                            degM::integer[4],
                            M::ARRAY(0..mLen-1, datatype=integer[8]),
                            uLen::integer[4],
                            degU::integer[4],
                            U::ARRAY(0..uLen-1, datatype=integer[8]),
                            N::integer[4],
                            DBound::integer[4],
                            p::integer[8],
                            nOLEN::integer[4],
                            nOUT::ARRAY(0..nOLEN-1, datatype=integer[8]),
                            degNOUT::REF(integer[4]),
                            dOLEN::integer[4],
                            dOUT::ARRAY(0..dOLEN-1, datatype=integer[8]),
                            degDOUT::REF(integer[4]),
                            RETURN::integer[4],
                            LIB=libR
);

# mRATRECON := subsop(1=(mLen,degM,M,uLen,degU,U,N,DBound,p,nOLEN,nOUT,degNOUT,dOLEN,dOUT,degDOUT),op(mRATRECON));

TIM := table():
CALLIDX := table():
CSVFILE := "calls.csv":
SUMMARYFILE := "summary.txt":

TimingReset := proc() option inline;
global TIM, CALLIDX, CSVFILE;
local fd;

    TIM["rr_cpp_total"] := 0.0:         TIM["rr_cpp_calls"] := 0:
    TIM["rr_cpp_ext_total"] := 0.0:     TIM["rr_cpp_ext_calls"] := 0:

    TIM["rr_map_total"] := 0.0:         TIM["rr_map_calls"] := 0:

    TIM["factor_total"] := 0.0:         TIM["factor_calls"] := 0:

    TIM["vsolve_cpp_total"] := 0.0:     TIM["vsolve_cpp_calls"] := 0:
    TIM["vsolve_cpp_ext_total"] := 0.0: TIM["vsolve_cpp_ext_calls"] := 0:

    TIM["vsolve_map_total"] := 0.0:     TIM["vsolve_map_calls"] := 0:

    CALLIDX["rr_cpp"] := 0:
    CALLIDX["rr_cpp_ext"] := 0:
    CALLIDX["rr_map"] := 0:
    CALLIDX["factor"] := 0:
    CALLIDX["vsolve_cpp"] := 0:
    CALLIDX["vsolve_cpp_ext"] := 0:
    CALLIDX["vsolve_map"] := 0:

    fd := fopen(CSVFILE, WRITE, TEXT):
    fprintf(fd, "event,call_index,j,tries,side,extra,dt\n"):
    fclose(fd):
end proc:

TIC := proc() option inline;
    return kernelopts(cputime);
end proc:

TADD := proc(keyTotal::string, keyCalls::string, dt::numeric) option inline;
global TIM;
    TIM[keyTotal] := TIM[keyTotal] + dt:
    TIM[keyCalls] := TIM[keyCalls] + 1:
end proc:

AVGSAFE := proc(total::numeric, calls::nonnegint) option inline;
    if calls = 0 then
        return 0.0;
    else
        return total/calls;
    fi;
end proc:

LogCSV := proc(event::string, j::{integer,string}, tries::{integer,string},
               side::string, extra::string, dt::numeric) option inline;
global CALLIDX, CSVFILE;
local fd, idx;

    CALLIDX[event] := CALLIDX[event] + 1:
    idx := CALLIDX[event]:

    fd := fopen(CSVFILE, APPEND, TEXT):
    fprintf(fd, "%s,%d,%a,%a,%s,%s,%g\n",
            event, idx, j, tries, side, extra, dt):
    fclose(fd):
end proc:

WriteTimings := proc(fname::string,
                     p::prime, T::posint, termsN::posint, termsD::posint,
                     N::nonnegint, DD::nonnegint,
                     numCount::posint, denCount::posint,
                     sameRF::{truefalse,boolean},
                     sameNum::{truefalse,boolean},
                     sameDen::{truefalse,boolean}) option inline;
global TIM, CSVFILE;
local fd;

    fd := fopen(fname,WRITE,TEXT):

    fprintf(fd,"====================================================\n"):
    fprintf(fd,"Timing summary in CPU seconds.\n"):
    fprintf(fd,"====================================================\n\n"):

    fprintf(fd,"Prime P                   -> = %a\n",p):
    fprintf(fd,"T (Assuming we know this) -> = %a\n",T):
    fprintf(fd,"Terms in Numerator        -> = %a\n",termsN):
    fprintf(fd,"Terms in Denominator      -> = %a\n",termsD):
    fprintf(fd,"Deg. numerator bound N    -> = %a\n",N):
    fprintf(fd,"Deg. denominator bound D  -> = %a\n",DD):
    fprintf(fd,"Numerator Count           -> = %a\n",numCount):
    fprintf(fd,"Denominator Count         -> = %a\n\n",denCount):

    fprintf(fd, "1) C++ wrapper RatRecon (Ratrecon1 total)\n"):
    fprintf(fd, "   calls = %d\n", TIM["rr_cpp_calls"]):
    fprintf(fd, "   total = %g\n", TIM["rr_cpp_total"]):
    fprintf(fd, "   avg   = %g\n\n",
            AVGSAFE(TIM["rr_cpp_total"], TIM["rr_cpp_calls"])):

    fprintf(fd, "   C++ external call only (mRATRECON)\n"):
    fprintf(fd, "   calls = %d\n", TIM["rr_cpp_ext_calls"]):
    fprintf(fd, "   total = %g\n", TIM["rr_cpp_ext_total"]):
    fprintf(fd, "   avg   = %g\n\n",
            AVGSAFE(TIM["rr_cpp_ext_total"], TIM["rr_cpp_ext_calls"])):

    fprintf(fd, "2) Maple RatRecon\n"):
    fprintf(fd, "   calls = %d\n", TIM["rr_map_calls"]):
    fprintf(fd, "   total = %g\n", TIM["rr_map_total"]):
    fprintf(fd, "   avg   = %g\n\n",
            AVGSAFE(TIM["rr_map_total"], TIM["rr_map_calls"])):

    fprintf(fd, "3) Maple Factor(annP) mod p\n"):
    fprintf(fd, "   calls = %d\n", TIM["factor_calls"]):
    fprintf(fd, "   total = %g\n", TIM["factor_total"]):
    fprintf(fd, "   avg   = %g\n\n",
            AVGSAFE(TIM["factor_total"], TIM["factor_calls"])):

    fprintf(fd, "4) C++ wrapper vSolve (VandermondeSolve1 total)\n"):
    fprintf(fd, "   calls = %d\n", TIM["vsolve_cpp_calls"]):
    fprintf(fd, "   total = %g\n", TIM["vsolve_cpp_total"]):
    fprintf(fd, "   avg   = %g\n\n",
            AVGSAFE(TIM["vsolve_cpp_total"], TIM["vsolve_cpp_calls"])):

    fprintf(fd, "   C++ external call only (mVSOLVE)\n"):
    fprintf(fd, "   calls = %d\n", TIM["vsolve_cpp_ext_calls"]):
    fprintf(fd, "   total = %g\n", TIM["vsolve_cpp_ext_total"]):
    fprintf(fd, "   avg   = %g\n\n",
            AVGSAFE(TIM["vsolve_cpp_ext_total"], TIM["vsolve_cpp_ext_calls"])):

    fprintf(fd, "5) Maple vSolve (VSolveMap)\n"):
    fprintf(fd, "   calls = %d\n", TIM["vsolve_map_calls"]):
    fprintf(fd, "   total = %g\n", TIM["vsolve_map_total"]):
    fprintf(fd, "   avg   = %g\n\n",
            AVGSAFE(TIM["vsolve_map_total"], TIM["vsolve_map_calls"])):

    fprintf(fd, "CSV per-call log file = %s\n\n", CSVFILE):

    fprintf(fd, "====================================================\n"):
    fprintf(fd, "Final correctness checks\n"):
    fprintf(fd, "====================================================\n"):
    fprintf(fd, "Same rational function?                    -> %a\n", sameRF):
    fprintf(fd, "Recovered numerator matches normalized?   -> %a\n", sameNum):
    fprintf(fd, "Recovered denominator matches normalized? -> %a\n", sameDen):

    fclose(fd):
    printf("Wrote timing summary to: %s\n", fname):
    printf("Wrote per-call CSV log to: %s\n", CSVFILE):
end proc:

TimingReset():

(* Michaels Code *)
VSolveMap := proc(m, v, p, shift::integer:=0 ) option inline;
local t,i,j,M,x,a,q,r,s;

   t := numelems(v);
   if numelems(m) <> t then
       error "v and m must be the same size";
   fi;

   printf("Maple code Vandermonde solver: t=%d  p=%d\n",t,p);

   M := 1;
   for r in m do
       M := Expand( M*(x-r) ) mod p;
   od;

   a := Vector(t);

   for j to t do
       q := Quo(M,x-m[j],x) mod p;
       r := 1/Eval(q,x=m[j]) mod p;
       s := 0;
       for i to t do
           s := s + v[i]*coeff(q,x,i-1);
       od;
       a[j] := r*s mod p;
       if shift <> 0 then
           r := 1/m[j] mod p;
           r := r &^ shift mod p;
           a[j] := r*a[j] mod p;
       fi;
   od;

   if type(v,list) then
       a := convert(a,list);
   fi;

   return a;
end proc:

VandermondeSolve1 := proc(m, v, p, shift::integer:=0) option inline;
local t,i,R,y,a,M,t0_ext,dt_ext;

    t := numelems(v);
    if numelems(m) <> t then
        error "v and m must be the same size";
    fi;

    R := Array(0..t-1,datatype=integer[8]);
    y := Array(0..t-1,datatype=integer[8]);
    a := Array(0..t-1,datatype=integer[8]);
    M := Array(0..t,datatype=integer[8]);

    for i from 0 to t-1 do
        R[i] := m[i+1];
        y[i] := v[i+1];
        a[i] := 0;
    od:

    for i from 0 to t do
        M[i] := 0;
    od:

    t0_ext := TIC():
    mVSOLVE(R,y,t,a,M,shift,p):
    dt_ext := TIC()-t0_ext:
    TADD("vsolve_cpp_ext_total","vsolve_cpp_ext_calls", dt_ext):
    LogCSV("vsolve_cpp_ext", "-", "-", "proc", sprintf("t=%d", t), dt_ext):

    return [seq(a[i], i=0..t-1)];
end proc:

POLYTOARR0 := proc(poly, var, deg, p) option inline;
local A,i;

    A := Array(0..deg, datatype=integer[8]);
    for i from 0 to deg do
        A[i] := coeff(poly,var,i) mod p;
    od:

    return A;
end proc:

ARRTOPOLYNEW := proc(A,deg,var) option inline;
local i:
    add(A[i]*var^i,i=0..deg);
end proc:

ARRTOPOLY0 := proc(A, deg, var, p)
local i, out;

    if deg < 0 then
        return 0;
    fi;

    out := 0;
    for i from 0 to deg do
        out := out + (A[i] mod p)*var^i;
    od:

    return expand(out) mod p;
end proc:

ArrayDegree0 := proc(A, len) option inline;
local i;
    for i from len-1 by -1 to 0 do
        if A[i] <> 0 then
            return i;
        fi;
    od;
    return -1;
end proc:

lastRatReconRC := 0:

Ratrecon1 := proc(Uin, Min, var,
                  N, DBound, p) option inline;
global lastRatReconRC;
local Upoly, Mpoly, degU, degM, uLen, mLen,
      UArr, MArr, nOLEN, dOLEN, nOUT, dOUT,
      degNOUT, degDOUT, rc, nn, dd, i, t0_ext, dt_ext;

    Upoly := Uin:
    Mpoly := Min:

    degU := degree(Upoly,var):
    degM := degree(Mpoly,var):

    if degU < 0 or degM < 0 then
        lastRatReconRC := -999:
        return FAIL;
    fi;

    uLen := degU+1:
    mLen := degM+1:

    UArr := POLYTOARR0(Upoly,var,degU,p):
    MArr := POLYTOARR0(Mpoly,var,degM,p):

    nOLEN := N+1:
    dOLEN := DBound+1:

    nOUT := Array(0..nOLEN-1, datatype=integer[8]):
    dOUT := Array(0..dOLEN-1, datatype=integer[8]):

    degNOUT := -1:
    degDOUT := -1:

    t0_ext := TIC():
    to 10^6 do:
    rc := mRATRECON(
        mLen, degM, MArr,
        uLen, degU, UArr,
        N, DBound, p,
        nOLEN, nOUT, degNOUT,
        dOLEN, dOUT, degDOUT
    ):
    od:
    dt_ext := TIC()-t0_ext:
    printf("DEGNOUT = %d, DEGDOUT = %d, TIME -> %.9f\n",degNOUT,degDOUT,dt_ext/10^6);
    quit;
    TADD("rr_cpp_ext_total","rr_cpp_ext_calls", dt_ext):
    LogCSV("rr_cpp_ext", "-", "-", "proc",
           cat("degU=", convert(degU,string), ";degM=", convert(degM,string)),
           dt_ext):
    
    lastRatReconRC := rc:

    if rc <> 0 then
        return FAIL;
    fi;

    if degNOUT < 0 or degNOUT > N then
        degNOUT := ArrayDegree0(nOUT, nOLEN);
    fi:

    if degDOUT < 0 or degDOUT > DBound then
        degDOUT := ArrayDegree0(dOUT, dOLEN);
    fi:

    if degNOUT < 0 or degNOUT > N then
        return FAIL;
    fi:

    if degDOUT < 0 or degDOUT > DBound then
        return FAIL;
    fi:

    nn := ARRTOPOLYNEW(nOUT,degNOUT,var):
    dd := ARRTOPOLYNEW(dOUT,degDOUT,var):

    if dd = 0 then
        return FAIL;
    fi:

    return nn/dd;
end proc:

p := prevprime(2^63-1):
# RF := rand():
n := x[1]^5+randpoly([seq(x[i],i=1..5)],terms=250,degree=5) mod p:
d := x[1]^5+randpoly([seq(x[i],i=1..5)],terms=250,degree=5) mod p:

f := n/d:

(* Assuming we know this. *)
termsN := nops(n):
termsD := nops(d):

T := max(termsN,termsD):

(* Assuming we know this. *)
N := degree(n):
DD := degree(d):

test := gcd(n,d):
test;

BLACKBOXF := proc(paramA::list(integer), p)
local sigmaVal,result,i;

    sigmaVal := paramA:
    result := Eval(f,[seq(x[i]=sigmaVal[i],i=1..nops(indets(f)))]) mod p:
    return result:
end proc:

COMPUTEBETA := proc()
local betaVal,i,r;

    r := rand(0..p-1):
    betaVal := table():

    for i from 1 to nops(indets(f))-1 do
        betaVal[i] := r():
    od:

    betaVal := convert(betaVal,list):
    return betaVal:
end proc:

GETT := proc(j,betaList)
local sigmaVal,alphaVal,TVal,i,BLACKBOXT,r,a,seen;

    r := rand(1..p-1):
    TVal := table():

    sigmaVal := [seq(ithprime(i)^j,i=1..nops(indets(f)))]:

    seen := table():
    alphaVal := []:
    while nops(alphaVal) < N+DD+1 do
        a := r():
        if not assigned(seen[a]) then
            seen[a] := true:
            alphaVal := [op(alphaVal), a]:
        fi:
    od:

    BLACKBOXT := proc(paramA::integer)
    local alphaX,result;

        alphaX := paramA:
        result := BLACKBOXF(
            [alphaX,
             seq(betaList[i-1]*alphaX-(betaList[i-1]*sigmaVal[1])+sigmaVal[i],
                 i=2..nops(indets(f)))],
            p
        ):
        return result:
    end proc:

    for i from 1 to (N+DD+1) do
        TVal[i] := BLACKBOXT(alphaVal[i]):
    od:

    TVal := convert(TVal,list):
    return alphaVal, TVal:
end proc:

betaList := COMPUTEBETA():

INTERSTEP := proc(betaList)
global lastRatReconRC;
local alphaVal,TVal,interpVal,ratReconVal,mapRatRecon,M,rr,i,j,
      tries,maxTries,ok,t0,dt;

    print(lastRatReconRC);

    ratReconVal := table():
    mapRatRecon := table():
    maxTries := 10:

    for j from 1 to 2*T do
        ok := false:

        for tries from 1 to maxTries while not ok do
            alphaVal,TVal := GETT(j,betaList):
            interpVal := Interp(alphaVal,TVal,z) mod p:
            M := [seq(z-alphaVal[i],i=1..nops(alphaVal))]:
            M := Expand(convert(M,`*`)) mod p:

            t0 := TIC():
            rr := Ratrecon1(interpVal,M,z,N,DD,p):
            dt := TIC()-t0:
            TADD("rr_cpp_total","rr_cpp_calls", dt):
            LogCSV("rr_cpp", j, tries, "loop", "Ratrecon1_total", dt):

            t0 := TIC():
            mapRatRecon[j] := Ratrecon(interpVal,M,z,N,DD) mod p:
            dt := TIC()-t0:
            TADD("rr_map_total","rr_map_calls", dt):
            LogCSV("rr_map", j, tries, "loop", "Maple_Ratrecon", dt):

            if rr <> FAIL then
                ratReconVal[j] := rr:
                ok := true:
            fi;
        od:

        if not ok then
            error "CPP ratRecon wrapper failed at j=%1 after %2 tries; last rc=%3",
                  j, maxTries, lastRatReconRC;
        fi;
    od:

    ratReconVal := convert(ratReconVal,list):
    mapRatRecon := convert(mapRatRecon,list):

    print(nops(ratReconVal)):
    print(nops(mapRatRecon)):

    return ratReconVal;
end proc:

ratReconVal := INTERSTEP(betaList):

ratNumer := table():
ratDenum := table():

for i from 1 to nops(ratReconVal) do
    ratNumer[i] := numer(ratReconVal[i]):
    ratDenum[i] := denom(ratReconVal[i]):
od:

ratNumer := convert(ratNumer,list):
ratDenum := convert(ratDenum,list):

ROOTTOMON := proc(r::posint)
local mon, i, pr, e, q, rr;

    rr := r:
    mon := 1:

    for i from 1 to nops(indets(f)) do
        pr := ithprime(i):
        e := 0:
        while irem(rr,pr,'q') = 0 do
            rr := q:
            e := e+1:
        od:
        mon := mon*x[i]^e:
    od:

    return mon:
end proc:

ROOTTOEXPS := proc(r::posint)
local exps, i, pr, e, q, rr;

    rr := r:
    exps := Array(1..nops(indets(f))):

    for i from 1 to nops(indets(f)) do
        pr := ithprime(i):
        e := 0:
        while irem(rr,pr,'q') = 0 do
            rr := q:
            e := e+1:
        od:
        exps[i] := e:
    od:

    return [seq(exps[i], i=1..nops(indets(f)))];
end proc:

COEFFFROMROOT := proc(poly, r::posint)
local c, exps, i;

    exps := ROOTTOEXPS(r):
    c := poly:

    for i from 1 to nops(exps) do
        c := coeff(c, x[i], exps[i]):
    od:

    return expand(c) mod p;
end proc:

BUILDSPARSE := proc(coeffs::list, roots::list)
local i, out;

    out := 0:
    for i from 1 to nops(coeffs) do
        out := out + coeffs[i]*ROOTTOMON(roots[i]):
    od:
    return expand(out) mod p:
end proc:

GETV := proc(paramA::list(polynom))
local listEval,result,j;

    result := table():
    listEval := paramA:

    for j from 1 to 2*T do
        result[j] := Eval(listEval[j],z=2^j) mod p:
    od:

    result := convert(result,list):
    return result:
end proc:

vNumer := GETV(ratNumer):
vDenom := GETV(ratDenum):

CREATEHANKELMAT := proc(paramA::list(integer), tdim::posint)
local val,i,j,H;

    val := paramA:
    H := Matrix(tdim,tdim):

    for i from 1 to tdim do
        for j from 1 to tdim do
            H[i,j] := val[i+j-1]:
        od:
    od:

    return H:
end proc:

GETBVEC := proc(paramA::list(integer), tdim::posint)
local values,bVec,i;

    values := paramA:
    bVec := -Vector([seq(values[i],i=tdim+1..2*tdim)]):
    return bVec:
end proc:

GETMONFACTOR := proc(paramA::Matrix, paramB::Vector, tdim::posint)
local HMat,bVec,L,annP,rootP,R,F,i,t0Factor,dtFactor;

    HMat := paramA:
    bVec := paramB:
    R := table():
    F := table():

    L := LinearSolve(HMat,bVec) mod p:

    annP := [z^tdim,seq(z^(i-1)*L[i],i=1..tdim)]:
    annP := add(annP[i],i=1..nops(annP)):

    t0Factor := TIC():
    annP := Factor(annP) mod p:
    dtFactor := TIC()-t0Factor:
    TADD("factor_total","factor_calls", dtFactor):
    LogCSV("factor", "-", "-", "proc", sprintf("tdim=%d", tdim), dtFactor):

    rootP := Roots(annP) mod p:

    for i from 1 to nops(rootP) do
        R[i] := rootP[i][1]:
    od:
    R := convert(R,list):

    for i from 1 to nops(R) do
        F[i] := ifactor(R[i]):
    od:
    F := convert(F,list):

    return R,F:
end proc:

HNumer := CREATEHANKELMAT(vNumer, termsN):
HDenom := CREATEHANKELMAT(vDenom, termsD):

bVecNumer := GETBVEC(vNumer, termsN):
bVecDenom := GETBVEC(vDenom, termsD):

numRoots,monNumFactor := GETMONFACTOR(HNumer,bVecNumer,termsN):
denomRoots,monDenomFactor := GETMONFACTOR(HDenom,bVecDenom,termsD):

numRoots := sort(numRoots):
denomRoots := sort(denomRoots):

numCount := nops(numRoots):
denCount := nops(denomRoots):

CREATEVANMAT := proc(paramA::list,n)
local xVals,i,j,vanMat;

    xVals := paramA:
    vanMat := Matrix(n,n):

    for i from 1 to n do
        for j from 1 to n do
            vanMat[i,j] := (xVals[j]^(i-1)) mod p:
        od:
    od:

    return vanMat:
end proc:

vanMatNumer := CREATEVANMAT(numRoots,numCount):
vanMatDenom := CREATEVANMAT(denomRoots,denCount):

GETBVECVAN := proc(paramA::list(integer),n)
local values,bVec,i;

    values := paramA:
    bVec := Vector([seq(values[i],i=1..n)]):
    return bVec:
end proc:

bVecVanNum := GETBVECVAN(vNumer,numCount):
bVecVanDenom := GETBVECVAN(vDenom,denCount):

NORMALIZEPAIR := proc(numCoeff::list, denCoeff::list, denRoots::list, p::prime)
local i, idx, lam, lamInv, newNum, newDen;

    idx := 0:

    for i from 1 to nops(denRoots) do
        if denRoots[i] = 1 then
            idx := i:
            break:
        fi:
    od:

    if idx = 0 then
        for i from 1 to nops(denCoeff) do
            if denCoeff[i] mod p <> 0 then
                idx := i:
                break:
            fi:
        od:
    fi:

    if idx = 0 then
        error "Cannot normalize -> all denominator coefficients are 0";
    fi:

    lam := denCoeff[idx] mod p:
    if lam = 0 then
        error "Cannot normalize -> all chosen denominator coefficient is 0";
    fi:

    lamInv := 1/lam mod p:

    newNum := [seq(lamInv*numCoeff[i] mod p, i=1..nops(numCoeff))]:
    newDen := [seq(lamInv*denCoeff[i] mod p, i=1..nops(denCoeff))]:

    return newNum, newDen, idx;
end proc:

NORMALIZEORIGINALPAIR := proc(nPoly, dPoly, denRoots::list, idx::posint, p::prime)
local lam, lamInv;

    lam := COEFFFROMROOT(dPoly, denRoots[idx]) mod p:
    if lam = 0 then
        error "Original denominator normalization coefficient is 0";
    fi:

    lamInv := 1/lam mod p:

    return expand(lamInv*nPoly) mod p, expand(lamInv*dPoly) mod p;
end proc:

t0 := TIC():
LVanNum := VandermondeSolve1( [seq(vNumer[i],i=1..numCount)],numRoots, p, 1):
dt := TIC()-t0:
TADD("vsolve_cpp_total","vsolve_cpp_calls", dt):
LogCSV("vsolve_cpp", "-", "-", "numer", sprintf("count=%d", numCount), dt):

t0 := TIC():
LVanDenom := VandermondeSolve1([seq(vDenom[i],i=1..denCount)],denomRoots, p, 1):
dt := TIC()-t0:
TADD("vsolve_cpp_total","vsolve_cpp_calls", dt):
LogCSV("vsolve_cpp", "-", "-", "denom", sprintf("count=%d", denCount), dt):

t0 := TIC():
LVanNumMap := VSolveMap([seq(vNumer[i],i=1..numCount)],numRoots, p, 1):
dt := TIC()-t0:
TADD("vsolve_map_total","vsolve_map_calls", dt):
LogCSV("vsolve_map", "-", "-", "numer", sprintf("count=%d", numCount), dt):

t0 := TIC():
LVanDenomMap := VSolveMap([seq(vDenom[i],i=1..denCount)],denomRoots, p, 1):
dt := TIC()-t0:
TADD("vsolve_map_total","vsolve_map_calls", dt):
LogCSV("vsolve_map", "-", "-", "denom", sprintf("count=%d", denCount), dt):

LVanNum, LVanDenom, normIdx := NORMALIZEPAIR(LVanNum, LVanDenom, denomRoots, p):

recNum := BUILDSPARSE(LVanNum, numRoots):
recDen := BUILDSPARSE(LVanDenom, denomRoots):

nNorm, dNorm := NORMALIZEORIGINALPAIR(n, d, denomRoots, normIdx, p):

sameRF := evalb(expand(recNum*d - recDen*n) mod p = 0):
sameNum := evalb(expand(recNum - nNorm) mod p = 0):
sameDen := evalb(expand(recDen - dNorm) mod p = 0):

printf("Same rational function? -> %a\n", sameRF):
printf("Recovered numerator matches normalized original? -> %a\n", sameNum):
printf("Recovered denominator matches normalized original? -> %a\n", sameDen):

WriteTimings(SUMMARYFILE,
             p, T, termsN, termsD, N, DD,
             numCount, denCount,
             sameRF, sameNum, sameDen):
