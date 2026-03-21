with(LinearAlgebra):

(* Returns a black box for the system of equations. *)
getBB := proc(arg1,arg2,arg3)
    local BB,L,sys,xVar,yVar:
    
    sys  := arg1:
    xVar := arg2:
    yVar := arg3:
    (*
    'L' is the coefficient matrix | b vector. In our case 
    the coefficients are in Z[x].
    *)
    L    := GenerateMatrix(sys,xVar,augmented=true):

    BB := proc(arg11::list(integer),p::prime)
        local pt,subVal,evalL,numEq,matA,errT,sol:

        pt := arg11:
        global counter:
        counter=counter+1:
        (*
        'subVal' is mapping each y[i] -> pt[i]. 
        *)
        subVal := zip((par,pnt) -> par = pnt,yVar,pt):
        evalL  := Eval(L,subVal) mod p:
        numEq  := numelems(xVar):
        (* Fill 'matA' with entries from 'evalL'. *)
        matA   := Matrix(numEq,numEq+1,datatype=integer[8],evalL):
        errT   := traperror(LinearSolve(p,matA,1)):
        if errT = "matrix is singular" then
            return FAIL:
        fi:
        sol := convert(matA[1..numEq,numEq+1],list):
        (* Solution vector. *)
        return sol:
    end proc:
    return BB:
end proc:



