getU := proc(M,col,alphaVal,p)
    local F,U,i:
    
    F := [seq(convert(M[..,i],list).i=1..col)]:
    U := [seq(Interp(alphaVal,F[i],x) mod p,i=1..col)]:
    return U:
end proc: