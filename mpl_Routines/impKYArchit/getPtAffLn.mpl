GPAFL := proc(numVar,alphaVal,betaVal,sigmaVal,p,T)
    option remember:
    local psi,nV,nP,i:
    
    # global numLines:
    # numLines := numLines+1:

    if numVar<2 then
        error "-1":
    fi:
    if nops(alphaVal)<T then
        error "-2":
    fi:
    for nP from 1 to T do
        psi[nP][1] := alphaVal[nP]:
        for nV from 2 to numVar do
            psi[nP][nV] := betaVal[nV-1]*alphaVal[nP]-betaVal[nV-1]*sigmaVal[1]+sigmaVal[nV] mod p:
        od:
    od:
    return [seq(convert(psi[i],list),i=1..T)]:
end proc: