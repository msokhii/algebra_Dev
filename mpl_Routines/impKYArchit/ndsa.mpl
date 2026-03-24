with(LinearAlgebra):

NDSA := proc(B,sigma_,beta_,num_var,p,num_points,num_eqn)
    local correct_degree,T,alpha,m,Psi_alpha,Y,u,dq,
          lin_sys,temp,result,count,M,row,col,DQ,MQRFR_done,i,r;

    print("In NDSA");

    MQRFR_done := Array(1..num_eqn, fill=false);
    result := Array(1..num_eqn);
    for i from 1 to num_eqn do
        result[i] := [];
    end do;

    correct_degree := false;
    lin_sys := false;
    T := num_points;
    count := 0;

    while true do
        count := count+1;
        print(count):
        print("T:= ",T);
        r := rand(p);
        alpha := [seq(i mod p, i=1..T)];
        lprint("NDSA: alpha: ",alpha);

        m := expand(product(x-alpha[j], j=1..T)) mod p;
        lprint("NDSA:m: ",m);

        Psi_alpha := GPAFL(num_var,alpha,beta_,sigma_,p,T);
        lprint("NDSA:Psi_alpha: ",Psi_alpha);

        Y := [seq(B(Psi_alpha[i],p), i=1..T)];
        print("NDSA: Y = ",Y);

        M := Matrix(Y);
        row, col := Dimension(M);
        lprint("row: ",row, " col: ",col);

        if row = 1 then
            lin_sys := false;
            u := Interp(alpha,Y,x) mod p;

            temp := [MQ(m,u,0,1,p)];
            result[1] := temp;
            dq := result[1][3];

        else
            lin_sys := true;
            u := getU(M,col,alpha,p);
        end if;

        if lin_sys = true then
            for i from 1 to nops(u) do
                if MQRFR_done[i] then
                    print("NDSA: Skipping MQRFR for equation ",i," as already done");
                    next;
                end if;

                temp := [MQ(m,u[i],0,1,p)];
                result[i] := temp;
                print("NDSA: result",i,": ",result[i]);
            end do;

            DQ := [seq(result[i][3], i=1..nops(u))];

            for i from 1 to nops(DQ) do
                if DQ[i] > 1 then
                    MQRFR_done[i] := true;
                    print("NDSA: MQRFR successful for equation ",i);
                end if;
            end do;

            dq := min(DQ);
        end if;

        if dq > 1 then
            print("NDSA: Termination condition met");
            return result, lin_sys;
        else
            print("NDSA:MQRFR failed. Trying again with more points");
            T := T*2;
            print("NDSA: mqrfr_status: ",MQRFR_done);

            for i from 1 to num_eqn do
                print("1"):
                if not MQRFR_done[i] then
                    print("2"):
                    result[i] := [];
                    print("3"):
                end if;
            end do;

            DQ := []; 
            print("4"):
        end if;

        if T > 2^5 then
            break;
        end if;
    end do;
end proc: