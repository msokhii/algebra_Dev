restart:
(* Helper Functions: *) 

(* ORDER: argA<-n,argB<-d,argC<-degN,argD<-degD *)
(* Function makes the denominator monic and scales the numerator 
   by the inverse of the denominator. *)

libNewton := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/newton.so":
libR := "/localhome/mss59/Desktop/research_Works/development/pol_ALGO/rfr.so":

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
                            ):

(* Remove type checking from mRATRECON. *)

mRATRECON := subsop(1=(
                       mLen,
                       degM,
                       M,
                       uLen,
                       degU,
                       U,
                       N,
                       DBound,
                       p,
                       nOLEN,
                       nOUT,
                       degNOUT,
                       dOLEN,
                       dOUT,
                       degDOUT),
                       op(mRATRECON)):

(* Converts a polynomial in Maple rep. to an array of coeffs. *)
convertPY2ARR := proc(poly,var,deg,p) option inline:
    local A,i;
    
    A := Array(0..deg,datatype=integer[8]):
    for i from 0 to deg do
        A[i] := coeff(poly,var,i):
    od:
    return A:
end proc:

(* Converts an array of coeffs. to a polynomial in Maple rep. 
   Array coeffs. are from low deg to high. *)
convertARR2PY := proc(A,deg,var) option inline:
    local i:
    
    add(A[i]*var^i,i=0..deg):
end proc:

(* Check if deg(A[i])=0 *)
checkZeroPY := proc(A,len) option inline:
    local i:

    for i from len-1 by -1 to 0 do
        if A[i]<>0 then
            return i:
        fi:
    od:
    printf("ZERO POLYNOMIAL.\n"):
    return -1:
end proc:

lastOP := 0:

(* Maple wrapper for mRATRECON. We call this maple function. *)
cppRR := proc(Uin,
              Min, 
              var,
              N, 
              DBound,
              p) option inline:
    
    global lastOP:
    local Upoly,Mpoly,degU,degM,uLen,mLen,
          UArr,MArr,nOLEN,dOLEN,nOUT,dOUT,
          degNOUT,degDOUT,cppRet,nn,dd,i:

        Upoly := Uin: #U
        Mpoly := Min: #M
        degU := degree(Upoly,var):
        degM := degree(Mpoly,var):

        if degU<0 or degM<0 then
            lastOP := -999: #U or M is a 0 polynomial.
            print("U OR M is a 0 polynomial : ",degU,degM):
            return FAIL:
        fi:

        uLen := degU+1:
        mLen := degM+1:

        (* Prepare inputs for cpp routine. *)
        UArr := convertPY2ARR(Upoly,var,degU,p):
        MArr := convertPY2ARR(Mpoly,var,degM,p):

        (* This is what the cpp routine will output. *)
        nOLEN := N+1:
        dOLEN := DBound+1:
        nOUT := Array(0..nOLEN-1,datatype=integer[8]):
        dOUT := Array(0..dOLEN-1,datatype=integer[8]):

        (* Initial State. Assuming N and D output will not be the 0 polynomial. *)
        degNOUT := -1:
        degDOUT := -1:
        
        cStart2 := time(): 
        to 10^3 do:
            cppRet := mRATRECON(
                                mLen,degM,MArr,
                                uLen,degU,UArr,
                                N,DBound,p,
                                nOLEN,nOUT,degNOUT,
                                dOLEN,dOUT,degDOUT
            ):
        od: 
        cStop2 := time()-cStart2: 
        printf("Local newton routine timing: %.9f\n",cStop2/10^3):
        lastOP := cppRet:
        
        (* 0 flag means success in reconstruction. *)
        if cppRet <> 0 then
            return FAIL:
        else
            print("RECONSTRUCTION SUCCESFULL!"):
        fi:
        if degNOUT<0 or degNOUT>N then
            degNOUT := checkZeroPY(nOUT,nOLEN);
        fi:
        if degDOUT<0 or degDOUT>DBound then
            degDOUT := checkZeroPY(dOUT,dOLEN);
        fi:
        if degNOUT<0 or degNOUT>N then
            print("degNOUT<0 or degNOUT>N!"):
            return FAIL:
        fi:
        if degDOUT<0 or degDOUT>DBound then
            print("degDOUT<0 or degDOUT>DBound!"):
            return FAIL:
        fi:
        
        (* Prepare output for maple rep. *)
        nn := convertARR2PY(nOUT,degNOUT,var):
        dd := convertARR2PY(dOUT,degDOUT,var):

        if dd=0 then
            print("DD=0 : ",dd):
            return FAIL:
        fi:

        (* Return reconstruction. *)
        return (nn/dd):
end proc:

mNEWTONINTERP := define_external(
                                'cppInterp',
                                xLen::integer[4],
                                xIn::ARRAY(0..xLen-1,datatype=integer[8]),
                                yLen::integer[4],
                                yIn::ARRAY(0..yLen-1,datatype=integer[8]),
                                p::integer[8],
                                outLen::integer[4],
                                yOut::ARRAY(0..outLen-1,datatype=integer[8]),
                                degOut::REF(integer[4]),
                                RETURN::integer[4],
                                LIB=libNewton
                                ):

(* Remove type checking from mNEWTONINTERP. *)

mNEWTONINTERP := subsop(1=(
                           xLen,
                           xIn,
                           yLen,
                           yIn,
                           p,
                           outLen,
                           yOut,
                           degOut),
                           op(mNEWTONINTERP)):

(* Check if deg(A[i])=0 *)
checkZeroPY := proc(A,len) option inline:
    local i:

    for i from len-1 by -1 to 0 do
        if A[i]<>0 then
            return i:
        fi:
    od:
    return -1:
end proc:

lastOP := 0:

cppNewtonInterp := proc(xVals,yVals,var,p) option inline:
    global lastOP;
    local n,xArr,yArr,outLen,yOut,degOut,cppRet,poly,i;
    
    (* Since input values are array's already. *)
    n := numelems(xVals):

    outLen := n:
    yOut := Array(0..outLen-1,datatype=integer[8]):
    degOut := 0:
    cStart := time():
    to 10^3 do
        cppRet := mNEWTONINTERP(
                                n,xVals,
                                n,yVals,
                                p,
                                outLen,yOut,degOut
        ):
    od: 
    cStop := time()-cStart: 
    printf("Local newton routine timing: %.9f\n",cStop/10^3):
    lastOP := cppRet:
    if cppRet <> 0 then
        return FAIL:
    else
        print("INTERPOLATION SUCCESFUL!"):
    fi:
    degOut := checkZeroPY(yOut,outLen):
    if degOut<0 then
        print("ZERO POLYNOMIAL : ",degOut):
        return FAIL:
    fi:
    poly := add(yOut[i]*var^i,i=0..degOut) mod p:
    return poly:
end proc:

MAKEMONIC := proc(argA::polynom,argB::polynom,argC::posint,argD::posint,p)
local n,d,g,degN,degD,LCD,invD:  

n := argA:
d := argB: 
LCD := lcoeff(d): 
if LCD<>1 then 
    g := Gcdex(LCD,p,x,'s','t') mod p:
    if g<>1 then 
        return 'FAIL':
    fi:
    invD := s:
    n := Expand(invD*n) mod p: 
    d := Expand(invD*d) mod p: 
    return (n,d): 
fi: 
return (n,d):
end proc:

(* ORDER: argA<-numPT *) 
(* Function generates values for x vector in preparation for 
   Newton Interpolation. *)
 
POPX := proc(argA::posint)
local xVec,numPT,i: 

numPT := argA:
xVec := Array(0..numPT-1,datatype=integer[8]):

for i from 0 to numPT-1 do 
    xVec[i] := i+1: 
od: 
return xVec:
end proc:     

EVALND := proc(argA::polynom,argB::polynom,argC::posint,p)
local tempDEval,tempNEval,invDEval,g,yArr,n,d,numPT,i: 

n := argA:
d := argB:
numPT := argC: 
yArr := Array(0..numPT-1,datatype=integer[8]): 

for i from 0 to numPT-1 do
    tempDEval := Eval(d,x=xArr[i]) mod p:
    if tempDEval=0 then
        return -1:
    fi:
    tempNEval := Eval(n,x=xArr[i]) mod p:
    g := Gcdex(tempDEval,p,x,'s','t') mod p:
    if g<>1 then
        return -2:
    fi: 
    invDEval := s: 
    yArr[i] := Expand(tempNEval*invDEval) mod p:
od: 
return yArr:
end proc:
(* Global Variables: *)

p := prevprime(2^32-1):
CT := 5:

(* Pseudo-random number generator: *)
prNum := rand(p):

(* Initial degrees: *)
degN := 5: 
degD := 5: 
for i from 1 to CT do
    printf("DEG N: %d DEG D: %d\n",degN,degD):
    n := x^degN+randpoly(x,coeffs=prNum,degree=degN) mod p:
    d := x^degD+randpoly(x,coeffs=prNum,degree=degD) mod p:
    n,d := MAKEMONIC(n,d,degN,degD,p): 
    f := n/d: 
    numPT := degN+degD+1:
    xArr := POPX(numPT):
    yArr := EVALND(n,d,numPT,p):
    tStart := time(): 
    to 10^3 do:
        mapU := Interp(xArr,yArr,z) mod p:
    od: 
    tStop := time()-tStart:
    printf("Maple newton routine timing: %.9f\n",tStop/10^3):
    localU := cppNewtonInterp(xArr,yArr,z,p):
    newtCheck := mapU-localU:
    print(newtCheck);
    M := [seq(z-xArr[i],i=1..numelems(xArr))]:
    M := Expand(convert(M,`*`)) mod p:
    tStart2 := time(): 
    to 10^3 do: 
        mapR := RatRecon(mapU,M,z,degN,degD) mod p: 
    od: 
    tStop2 := time()-tStart2: 
    printf("Maple ratrecon routine timing: %.9f\n",tStop2/10^3):
    localRR := cppRR(localU,M,z,degN,degD,p):
    rrCheck := mapR-localRR:
    print(rrCheck);
    degN := degN*2:
    degD := degD*2: 
od:
