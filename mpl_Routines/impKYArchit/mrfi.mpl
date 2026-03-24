MRFI := proc(BB,arg2,arg3,arg4,p)
    local numVars,numEq,vars,randNum,pList,sigmaVal,
          numEval,denomEval,T,T_old,ptMQRFR,numList,
          denumList,dir,linSys:

    numVars := arg2: (* Nops yVars. *)
    numEq   := arg3: (* Nops xVars. *)
    vars    := arg4: (* Actual yVars. *)

    lprint("STARTING MRFI: "):
    lprint("No. of parameters -> ",numVars):
    lprint("No. of equations  -> ",numEq):
    printf("\n\n"):

    randNum := rand(p):
    dir := [seq(randNum(),i=1..numVars-1)]:
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
    MQRFRres := NDSA(BB,[seq(1,i=1..numVars)],dir,numVars,p,T,numEq):
    lprint("IFDHSHFSDF"):
    print(MQRFRres):
    print(linSys):

    numList := [seq(MQRFRres[i][1],i=1..nops(MQRFRres))]:
    denumList := [seq(MQRFRres[i][2],i=1..nops(MQRFRres))]:
    print(numList);
    print(denumList);
    
    return numList,denumList:
end proc:
