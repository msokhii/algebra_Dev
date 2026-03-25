restart:

libR := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/rfr.so":

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

mRATRECON := subsop(1=(mLen,degM,M,uLen,degU,U,N,DBound,p,nOLEN,nOUT,degNOUT,dOLEN,dOUT,degDOUT),op(mRATRECON)):

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

cppRR := proc(Uin,
              Min, 
              var,
              N, 
              DBound,
              p) option inline;
    
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

        rc := mRATRECON(
            mLen, degM, MArr,
            uLen, degU, UArr,
            N, DBound, p,
            nOLEN, nOUT, degNOUT,
            dOLEN, dOUT, degDOUT
        ):
        
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

(* Inputs: *)

u := 3+5*x+7*x^2+11*x^3+13*x^4:
m := x^5:
p := 97:
N := 2:
DD := 2: 

mapRR := Ratrecon(u,m,x,N,DD) mod p:
print(mapRR):

myRR := cppRR(u,m,x,N,DD,p):
print(myRR):
print(lastRatReconRC):

