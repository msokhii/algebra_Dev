MRFI := proc(BB,arg2,arg3,arg4,p)
    local numVars,numEq,vars,randNum,pList,sigmaVal,
          numEval,denomEval,T,T_old,ptMQRFR,numList,
          denumList:

    numVars := arg2: (* Nops xVars. *)
    numEq   := arg3: (* Nops yVars. *)
    vars    := arg4: (* Actual yVars. *)

    lprint("STARTING MRFI: "):
    lprint("No. of parameters -> ",numVars):
    lprint("No. of equations  -> ",numEq):
    printf("\n\n"):

    randNum := rand(p):
    pList := [seq(ithprime(i),i=1..numVars)]:
    sigmaVal := []:
    
    (*
    List of lists. 
    *)
    numEval := [seq([],i=1..numEq)]:
    denEval := [seq([],i=1..numEq)]: 
    T       := 4:
    T_old   := 1:
    ptMQRFR := [seq(0,i=1..numEq)]:

    (* Initial NDSA call. *)
    MQRFRres,linSys := NDSA(BB,[seq(1,i=1..numVars)],[7],numVars,p,T,numEq):
    MQRFRres;
    linSys;

    numList := [seq(MQRFRres[i][1],i=1..nops(MQRFRres))]:
    denumList := [seq(MQRFRres[i][2],i=1..nops(MQRFRres))]:
    print(numList);
    print(denumList);
    
    return numList,denumList:
end proc:
