
genSystem := proc(arg::string)
    local system,vars,params,i:

    if arg = "bspline" then
        (* This is a 21x21 system. *)
        system := {x[7]+x[12]-1, 
                  x[8]+x[13]-1, 
                  x[21]+x[6]+x[11]-1, 
                  x[1]*y[1]+x[1]-x[2],
                  x[11]*y[3]+x[11]-x[12], 
                  x[16]*y[5]-x[17]*y[5]-x[17], 
                  -x[20]*y[3]+x[21]*y[3]+x[21], 
                  x[3]*y[2]+x[3]-x[4],
                  -x[8]*y[4]+x[9]*y[3]+x[9], 
                  2*x[1]*y[1]^2-2*x[1]-2*x[10]+4*x[2], 
                  -x[10]*y[2]+x[18]*y[2]+x[18]-x[19], 
                  2*x[11]*y[3]^2-2*x[11]+4*x[12]-2*x[13], 
                  -x[13]*y[4]+x[14]*y[4]+x[14]-x[15], 
                  2*x[15]*y[5]^2-4*x[16]*y[5]^2+2*x[17]*y[5]^2-2*x[17],
                  2*x[19]*y[3]^2-4*x[20]*y[3]^2+2*x[21]*y[3]^2-2*x[21], 
                  2*x[3]*y[2]^2-2*x[3]+4*x[4]-2*x[5], 
                  -x[5]*y[3]+x[6]*y[3]+x[6]-x[7], 
                  2*x[7]*y[4]^2-4*x[8]*y[4]^2+2*x[9]*y[4]^2-2*x[9], 
                  -4*x[10]*y[2]^2+2*x[18]*y[2]^2+2*x[2]*y[2]^2-2*x[18]+4*x[19]-2*x[20],
                  2*x[12]*y[4]^2-4*x[13]*y[4]^2+2*x[14]*y[4]^2-2*x[14]+4*x[15]-2*x[16], 
                  2*x[4]*y[3]^2-4*x[5]*y[3]^2+2*x[6]*y[3]^2-2*x[6]+4*x[7]-2*x[8]}:
    fi:

    (* 
    'vars' is the x variable in the system. 
    'params' is the y variable in the system.
    *)
    vars := {seq(x[i],i=1..nops(system))}:
    (*
    'params' peforms a set difference i.e all the indets which in our case 
    is 26 - 21 (number of x vars.) and outputs the y vars which is 5.
    *)
    params := indets(system) minus vars:
    return (system,convert(vars,list),convert(params,list),nops(params),nops(vars)):
end proc: