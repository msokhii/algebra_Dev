MQ := proc(r0::polynom,r1::polynom,t0::integer,t1::integer,p::prime)
    local r,t,q,i,f,g,qMax,lcG; 
    r[0] := r0:
    r[1] := r1:
    t[0] := t0:
    t[1] := t1:
    f    := r0:
    g    := t1:
    qMax := 0:
    i    := 1:
    while r[i] <> 0 do
        q[i] := Quo(r[i-1],r[i],x,'r[i+1]') mod p:
        if degree(q[i],x)>qMax then
            qMax := degree(q[i],x):
            f    := r[i]:
            g    := t[i]:
        fi:
        t[i+1] := Expand(t[i-1]-q[i]*t[i]) mod p:
        i      := i+1:
        if qMax <= 1 or gcd(f,g) <> 1 or g = 0 then
            FAIL:
        fi:
    od:
    lcG := lcoeff(g):
    return ((f/lcG) mod p,(g/lcG) mod p,qMax,lcG):
end proc: