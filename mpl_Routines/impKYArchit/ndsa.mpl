with(LinearAlgebra):

NDSA := proc(BB,sigmaVal,betaVal,numVar,p,numPt,numEq)
    local MQRFRDone,correctDeg,linSys,T,res,count,r,
          alphaVal,m,psiAlpha,yVal,M,row,col,u,
          degQ,temp: 

    printf("IN NDSA NOW!"):
    MQRFRDone  := [seq(false,i=1..numEq)]:
    correctDeg := false:
    linSys     := false:
    T          := numPt:
    temp       := []:
    res        := [seq([],i=1..numEq)]:
    count      := 0:

    while(not(correctDeg)) do
        count := count+1:
        printf("T: ",T):
        r := rand(p):
        alphaVal := [seq(i mod p,i=1..T)]:
        print(alphaVal):
        m := Expand(product(x-alphaVal[j],j=1..T)) mod p:
        print(betaVal):
        psiAlpha := GPAFL(numVar,alphaVal,betaVal,sigmaVal,p,T):
        print(psiAlpha):
        yVal := [seq(BB(psiAlpha[i],p),i=1..T)]:
        M := Matrix(yVal):
        row,col := Dimension(M):

        if row=1 then
            linSys := false:
            u      := Interp(alphaVal,yVal,x) mod p:
            res    := [[MQ(m,u,0,1,p)]]:
            print(res);
            degQ   := res[1][3]:
        else:
            linSys := true:
            u      := getU(M,col,alphaVal,p):
        fi:

        if linSys = true then 
            for i from 1 to nops(u) do
                if(MQRFRDone[i]) then
                    print("SKIP!"):
                    next:
                fi:
                temp := [op(temp),MQ(m,u[i],0,1,p)]:
                res[i] := temp:
                print(res);
                temp := []:
            od:
            DQ := [seq(res[i][3],i=1..nops(res))]:
            for i from 1 to numelems(DQ) do
                if DQ[i]>1 then
                    MQRFRDone[i] := true:
                fi:
            od:
            degQ := min(DQ):
        fi: 

        if degQ>1 then
            return (res,linSys):
        else
            T := T*2:
            for i from 1 to numEq do
                if MQRFRDone[i] = false then
                    res[i] := []:
                    print(res);
                fi:
            od:
            DQ := []:
        fi:
        
        if(T>2^5) then
            break:
        fi:
    od:
end proc: