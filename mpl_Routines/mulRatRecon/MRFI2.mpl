(*
MRFI2 := proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, Primes, direction, sigma_, num_eval, den_eval, u, mon,
          numerator_done, denominator_done, T, T_old,
          mqrfr_results, lin_sys, num_points_mqrfr,
          Numerators, Denominators, deg_num, deg_den,
          lambda_num, lambda_den, terms_num, terms_den,
          R_num, R_den, Roots_num_eval, Roots_den_eval,
          num_mono, den_mono, coeff_num, coeff_den,
          final_num, final_den, r_, temp, common_den_flag,
          den_lc, den_lc_inv, bmea_done, r, temp_den,
          all_den_done, all_done, max_num_points_mqrfr,
          tempTable, ck, tempN, tempD, rTemp, alpha, Psi_alpha,
          Y, U, m, iter, iter2, mqrfr_list;

    lprint("MRFI2  ========================================"):
    lprint("MRFI2 Starting MRFI2"):
    lprint("MRFI2 Number of parameters:", num_vars):
    lprint("MRFI2 Number of equations:", num_eqn):
    lprint("MRFI2 ========================================"):

    r_ := rand(p):
    Primes := [seq(ithprime(i), i=1..num_vars)]:
    direction := [seq(r_(), i=1..num_vars-1)]:
    sigma_ := []:

    # Initializing accumulators and flags
    num_eval := [seq([], i=1..num_eqn)]:
    den_eval := [seq([], i=1..num_eqn)]:
    
    numerator_done := [seq(false, i=1..num_eqn)]:
    denominator_done := [seq(false, i=1..num_eqn)]:
    bmea_done := [seq(false, i=1..num_eqn)]:
    all_done := true:

    lambda_num := table():
    terms_num := table():
    R_num := table():
    lambda_den := table():
    terms_den := table():
    R_den := table():
    final_num := table():
    final_den := table():
    num_mono := table():
    coeff_num := table():
    den_mono := table():
    coeff_den := table():

    for i from 1 to num_eqn do
        lambda_num[i] := []:
        terms_num[i] := []:
        R_num[i] := []:
        lambda_den[i] := []:
        terms_den[i] := []:
        R_den[i] := []:
    end do:

    lprint("MRFI2 numerator_done: ", numerator_done):
    lprint("MRFI2 denominator_done: ", denominator_done):

    T := 4:
    T_old := 1:
    
    print(T):
    mqrfr_results,T,lin_sys:= NDSA(B, [seq(1, i=1..num_vars)], direction, num_vars, p, T, num_eqn):
    # print(mqrfr_results):
    print(T):

    Numerators := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    Denominators := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
    print(Numerators):
    print(Denominators):

    for k from 1 to num_eqn do
        num_eval[k] := [eval(Numerators[k], x=1) mod p]:
        den_eval[k] := [eval(Denominators[k], x=1) mod p]:
    end do:

    deg_num := [seq(degree(Numerators[i], x), i=1..nops(Numerators))]:
    deg_den := [seq(degree(Denominators[i], x), i=1..nops(Denominators))]:

    common_den_flag := true:
    for k from 2 to num_eqn do
        if Denominators[k] <> Denominators[1] then
            common_den_flag := false:
            break:
        end if:
    end do:
    # lprint("MRFI2 Common denominator:", common_den_flag):
    print("______________________________________________________________________________"):
    print("in main evaluation loop of MRFI2"):

    print(deg_num): 
    print(deg_den):  
    print(nops(deg_num)):
    print(nops(deg_den)):
    print(max(deg_num)):
    print(max(deg_den)):

    tempG := rand(p):
    ratReconVal := table():

    for i from 1 to nops(deg_num) do
        tempDegNum := deg_num[i]:
        tempDegDenom := deg_den[i]:
        ratReconVal[i] := table():

        for j from 1 to 2*T do
            sigma_j := [seq(Primes[k]^j mod p, k=1..nops(Primes))]:
            alphaVal := [seq(tempG(), s=1..T)]:
            m := [seq(x-alphaVal[s], s=1..nops(alphaVal))]:
            m := Expand(convert(m,`*`)) mod p:

            Psi_alpha := get_point_on_affine_line(num_vars, alphaVal, direction, sigma_j, p, T):
            Y := [seq(B(Psi_alpha[s],p)[i], s=1..nops(alphaVal))]:
            interpVal := Interp(alphaVal, Y, x) mod p:

            rr := Ratrecon(interpVal, m, x, tempDegNum, tempDegDenom) mod p:
            ratReconVal[i][j] := rr:

            num_eval[i] := [op(num_eval[i]), eval(numer(rr), x=sigma_j[1]) mod p]:
            den_eval[i] := [op(den_eval[i]), eval(denom(rr), x=sigma_j[1]) mod p]:
        od:
    od:

    numList := [seq([seq(numer(ratReconVal[i][j]), j=1..2*T)], i=1..nops(deg_num))]:
    denList := [seq([seq(denom(ratReconVal[i][j]), j=1..2*T)], i=1..nops(deg_num))]:
    print("MRFI2 num_eval: ",num_eval):
    print("MRFI2 den_eval: ",den_eval):
    print(nops(num_eval)):
    print(nops(den_eval)):
    print("______________________________________________________________________________"):
    print("numerator_done: ",numerator_done):
    print("denominator_done: ",denominator_done):
    print("--------------------------------------------------------------------------------"):

    all_done := true:
        for k from 1 to num_eqn do
            lprint("MRFI2 k=",k):
            lprint("MRFI2 numerator_done[",k,"]=",numerator_done[k]):
            if  numerator_done[k] then print("MRFI2 Skipping BMEA for numerator component ",k," as already done"): next: end if;
                # temp := BMEA(num_eval[k], p, Z):
                # lprint("MRFI2 temp (numerator)=",temp):
                
                # Use hash table
                lambda_num[k] := BMEA(num_eval[k], p, Z):
                terms_num[k] := degree(lambda_num[k], Z):
                R_num[k] := Roots(lambda_num[k]) mod p:

                lprint("MRFI2 lambda_num: ",lambda_num[k]):
                lprint("MRFI2 terms_num: ",terms_num[k]):
                lprint("MRFI2 R_num: ",R_num[k]):
                print("--------------------------------------------------------------------------------"):
                print("Checking termination condition for numerator"):
                print("nops(R_num[",k,"]): ", nops(R_num[k])):
                print("R_num[",k,"]: ", R_num[k]);
                print("terms_num[",k,"]: ", terms_num[k]);
                print("T: ", T);
                
                # Add check for empty R_num[k]
                if R_num[k] = [] then
                    print("MRFI2: Empty roots list for numerator component ", k):
                    numerator_done[k] := false:
                    next:
                end if:

                if nops(R_num[k]) > 0 and R_num[k][1][1] = 0 then
                    R_num[k] := remove(x->x=[0,1], R_num[k]):
                    terms_num[k] := terms_num[k] - 1:
                end if:

                if nops(R_num[k]) = terms_num[k] and terms_num[k] <= T then
                    print("MRFI2 Numerator component ",k," recovered!"):
                    numerator_done[k] := true:
                end if:
            print("____________________________________________________________________________"):
        end do:

        print("===================================================================================="):
        print("MFRI processing denominators now"):

        # Process denominators
        lprint("MRFI2 den_eval: ",den_eval):

        all_den_done := true:
        lprint("MRFI2 common_den_flag: ",common_den_flag):
        if common_den_flag then
            lprint("In common_den_flag"):
            if not denominator_done[1] then
                
                # temp := BMEA(den_eval[1], p, Z):
                # lprint("MRFI2 temp (denominator)=",temp):
                lambda_den[1] := BMEA(den_eval[1], p, Z):
                terms_den[1] :=degree(lambda_den[1], Z):
                R_den[1] := Roots(lambda_den[1]) mod p:
                lprint("MRFI2 lambda_den: ",lambda_den[1]):
                lprint("MRFI2 terms_den: ",terms_den[1]):
                lprint("MRFI2 R_den: ",R_den[1]):
                
                if nops(R_den[1]) > 0 and R_den[1][1][1] = 0 then
                    R_den[1] := remove(x->x=[0,1], R_den[1]):
                    terms_den[1] := terms_den[1] - 1:
                end if:
                
                if nops(R_den[1]) = terms_den[1] and terms_den[1] < T then
                    for k from 1 to num_eqn do
                        denominator_done[k] := true:
                    end do:
                end if:
            end if:
        else
            for k from 1 to num_eqn do
                if  denominator_done[k] then 
                    print("MRFI2 Skipping BMEA for denominator component ",k," as already done"): 
                    next:
                end if;
                # temp := BMEA(den_eval[k], p, Z):
                
                # Use hash table
                lambda_den[k] := BMEA(den_eval[k], p, Z):
                terms_den[k] := degree(lambda_den[k], Z):
                R_den[k] := Roots(lambda_den[k]) mod p:

                    print("--------------------------------------------------------------------------------"):
                    print("Checking termination condition for denominator"):
                    print("terms_den[",k,"]: ", terms_den[k]);
                    print("T: ", T);
                    
                    if nops(R_den[k]) > 0 and R_den[k][1][1] = 0 then
                        R_den[k] := remove(x->x=[0,1], R_den[k]):
                        terms_den[k] := terms_den[k] - 1:
                    end if:
                    
                    if nops(R_den[k]) = terms_den[k] and terms_den[k] < T then
                        print("MRFI2 Denominator component ",k," recovered!"):
                        denominator_done[k] := true:
                    end if:
            end do:
        end if:

        print("--------------------------------------------------------------------------------"):
        for i from 1 to num_eqn do 
            bmea_done[i] := numerator_done[i] and denominator_done[i]:
        end do:
        print("BMEA done status: ",bmea_done):
        for i from 1 to num_eqn do 
            all_done:=all_done and bmea_done[i]:
        end do:
        print("All done status: ",all_done):
        if all_done then
            lprint("MRFI2 All components recovered!"):
        end if:
        print("numerator_done: ",numerator_done):
        print("denominator_done: ",denominator_done):
        print("R_num: ",R_num):
        print("R_den: ",R_den):
        print("terms_num: ",terms_num):
        print("terms_den: ",lambda_num):
        print("lambda_den: ",lambda_den):

end proc:
*)

MRFI2 := proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, Primes, direction, sigma_j, num_eval, den_eval,
          numerator_done, denominator_done, T,T_old,
          mqrfr_results, lin_sys,
          Numerators, Denominators, deg_num, deg_den,
          lambda_num, lambda_den, terms_num, terms_den,
          R_num, R_den,
          final_num, final_den,
          coeff_num, coeff_den,
          common_den_flag,
          bmea_done, all_done,
          r_, tempG, alphaVal, Psi_alpha, Y, m, rr,
          ratReconVal, numList, denList,
          tempDegNum, tempDegDenom, sampleCount;

    lprint("MRFI2 ========================================"):
    lprint("MRFI2 Starting MRFI2"):
    lprint("MRFI2 Number of parameters:", num_vars):
    lprint("MRFI2 Number of equations:", num_eqn):
    lprint("MRFI2 ========================================"):

    r_ := rand(p):
    Primes := [seq(ithprime(i), i=1..num_vars)]:
    direction := [seq(r_(), i=1..num_vars-1)]:

    numerator_done := [seq(false, i=1..num_eqn)]:
    denominator_done := [seq(false, i=1..num_eqn)]:
    bmea_done := [seq(false, i=1..num_eqn)]:

    lambda_num := table():
    terms_num := table():
    R_num := table():
    lambda_den := table():
    terms_den := table():
    R_den := table():
    final_num := table():
    final_den := table():
    coeff_num := table():
    coeff_den := table():

    for i from 1 to num_eqn do
        lambda_num[i] := []:
        terms_num[i] := 0:
        R_num[i] := []:
        lambda_den[i] := []:
        terms_den[i] := 0:
        R_den[i] := []:
    end do:

    T := 4:
    mqrfr_results,T,lin_sys := NDSA(B, [seq(1, i=1..num_vars)], direction, num_vars, p, T, num_eqn):

    Numerators := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    Denominators := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:

    deg_num := [seq(degree(Numerators[i], x), i=1..nops(Numerators))]:
    deg_den := [seq(degree(Denominators[i], x), i=1..nops(Denominators))]:

    common_den_flag := true:
    for k from 2 to num_eqn do
        if Denominators[k] <> Denominators[1] then
            common_den_flag := false:
            break:
        end if:
    end do:

    print("______________________________________________________________________________"):
    print("in main evaluation loop of MRFI2"):
    print("T returned by NDSA: ", T):
    print("deg_num = ", deg_num):
    print("deg_den = ", deg_den):

    tempG := rand(p):
    ratReconVal := table():

    # Keep the initial sample at x = 1
    num_eval := [seq([eval(Numerators[k], x=1) mod p], k=1..num_eqn)]:
    den_eval := [seq([eval(Denominators[k], x=1) mod p], k=1..num_eqn)]:

    # Collect additional samples for j = 1 .. 2*T-1
    T := 32:
    for j from 1 to 2*T-1 do
        sigma_j := [seq(Primes[k]^j mod p, k=1..nops(Primes))]:

        for i from 1 to nops(deg_num) do
            tempDegNum := deg_num[i]:
            tempDegDenom := deg_den[i]:

            if not assigned(ratReconVal[i]) then
                ratReconVal[i] := table():
            end if:

            alphaVal := [seq(tempG(), s=1..T)]:
            m := Expand(mul(x-alphaVal[s], s=1..nops(alphaVal))) mod p:

            Psi_alpha := get_point_on_affine_line(num_vars, alphaVal, direction, sigma_j, p, T):
            Y := [seq(B(Psi_alpha[s], p)[i], s=1..nops(alphaVal))]:

            interpVal := Interp(alphaVal, Y, x) mod p:
            rr := Ratrecon(interpVal, m, x, tempDegNum, tempDegDenom) mod p:
            ratReconVal[i][j] := rr:

            num_eval[i] := [op(num_eval[i]), Eval(numer(rr), x=sigma_j[1]) mod p]:
            den_eval[i] := [op(den_eval[i]), Eval(denom(rr), x=sigma_j[1]) mod p]:
        end do:
    end do:

    numList := [seq([seq(numer(ratReconVal[i][j]), j=1..2*T-1)], i=1..nops(deg_num))]:
    denList := [seq([seq(denom(ratReconVal[i][j]), j=1..2*T-1)], i=1..nops(deg_num))]:

    print("MRFI2 num_eval: ", num_eval):
    print("MRFI2 den_eval: ", den_eval):
    print("______________________________________________________________________________"):
    print("numerator_done: ",numerator_done):
    print("denominator_done: ",denominator_done):
    print("--------------------------------------------------------------------------------"):

    all_done := true:
        for k from 1 to num_eqn do
            lprint("MRFI k=",k):
            lprint("MRFI numerator_done[",k,"]=",numerator_done[k]):
            if  numerator_done[k] then print("MRFI Skipping BMEA for numerator component ",k," as already done"): next: end if;
                # temp := BMEA(num_eval[k], p, Z):
                # lprint("MRFI temp (numerator)=",temp):
                
                # Use hash table
                lambda_num[k] := BMEA(num_eval[k], p, Z):
                terms_num[k] := degree(lambda_num[k], Z):
                R_num[k] := Roots(lambda_num[k]) mod p:

                lprint("MRFI lambda_num: ",lambda_num[k]):
                lprint("MRFI terms_num: ",terms_num[k]):
                lprint("MRFI R_num: ",R_num[k]):
                print("--------------------------------------------------------------------------------"):
                print("Checking termination condition for numerator"):
                print("nops(R_num[",k,"]): ", nops(R_num[k])):
                print("R_num[",k,"]: ", R_num[k]);
                print("terms_num[",k,"]: ", terms_num[k]);
                print("T: ", T);
                
                # Add check for empty R_num[k]
                if R_num[k] = [] then
                    print("MRFI: Empty roots list for numerator component ", k):
                    numerator_done[k] := false:
                    next:
                end if:

                if nops(R_num[k]) > 0 and R_num[k][1][1] = 0 then
                    R_num[k] := remove(x->x=[0,1], R_num[k]):
                    terms_num[k] := terms_num[k] - 1:
                end if:

                if nops(R_num[k]) = terms_num[k] and terms_num[k] <= T then
                    print("MRFI Numerator component ",k," recovered!"):
                    numerator_done[k] := true:
                end if:
            print("____________________________________________________________________________"):
        end do:
        print("===================================================================================="):
        print("MFRI processing denominators now"):
        
        # Process denominators
        lprint("MRFI den_eval: ",den_eval):
        
        all_den_done := true:
        lprint("MRFI common_den_flag: ",common_den_flag):
        if common_den_flag then
            lprint("In common_den_flag"):
            if not denominator_done[1] then
                
                # temp := BMEA(den_eval[1], p, Z):
                # lprint("MRFI temp (denominator)=",temp):
                lambda_den[1] := BMEA(den_eval[1], p, Z):
                terms_den[1] :=degree(lambda_den[1], Z):
                R_den[1] := Roots(lambda_den[1]) mod p:
                lprint("MRFI lambda_den: ",lambda_den[1]):
                lprint("MRFI terms_den: ",terms_den[1]):
                lprint("MRFI R_den: ",R_den[1]):
                
                if nops(R_den[1]) > 0 and R_den[1][1][1] = 0 then
                    R_den[1] := remove(x->x=[0,1], R_den[1]):
                    terms_den[1] := terms_den[1] - 1:
                end if:
                
                if nops(R_den[1]) = terms_den[1] and terms_den[1] < T then
                    for k from 1 to num_eqn do
                        denominator_done[k] := true:
                    end do:
                end if:
            end if:
        else
            for k from 1 to num_eqn do
                if  denominator_done[k] then 
                    print("MRFI Skipping BMEA for denominator component ",k," as already done"): 
                    next:
                end if;
                # temp := BMEA(den_eval[k], p, Z):
                
                # Use hash table
                lambda_den[k] := BMEA(den_eval[k], p, Z):
                terms_den[k] := degree(lambda_den[k], Z):
                R_den[k] := Roots(lambda_den[k]) mod p:

                    print("--------------------------------------------------------------------------------"):
                    print("Checking termination condition for denominator"):
                    print("terms_den[",k,"]: ", terms_den[k]);
                    print("T: ", T);
                    
                    if nops(R_den[k]) > 0 and R_den[k][1][1] = 0 then
                        R_den[k] := remove(x->x=[0,1], R_den[k]):
                        terms_den[k] := terms_den[k] - 1:
                    end if:
                    
                    if nops(R_den[k]) = terms_den[k] and terms_den[k] < T then
                        print("MRFI Denominator component ",k," recovered!"):
                        denominator_done[k] := true:
                    end if:
            end do:
        end if:
        print(denominator_done):
        print(numerator_done):

        print("--------------------------------------------------------------------------------"):
        for i from 1 to num_eqn do 
            bmea_done[i] := numerator_done[i] and denominator_done[i]:
        end do:
        print("BMEA done status: ",bmea_done):
        for i from 1 to num_eqn do 
            all_done:=all_done and bmea_done[i]:
        end do:
        print("All done status: ",all_done):
        if all_done then
            lprint("MRFI All components recovered!"):
        end if:
        print("numerator_done: ",numerator_done):
        print("denominator_done: ",denominator_done):
        print("R_num: ",R_num):
        print("R_den: ",R_den):
        print("terms_num: ",terms_num):
        print("terms_den: ",lambda_num):
        print("lambda_den: ",lambda_den):

        Roots_num_eval := [seq([seq(r[1], r in R_num[k])], k=1..num_eqn)]:
    print("MRFI Roots_num_eval: ",Roots_num_eval):
    
    if common_den_flag then
        Roots_den_eval := [[seq(r[1], r in R_den[1])]]:
        for k from 2 to num_eqn do
            Roots_den_eval := [op(Roots_den_eval), Roots_den_eval[1]]:
        end do:
    else
        Roots_den_eval := [seq([seq(r[1], r in R_den[k])], k=1..num_eqn)]:
    end if:
    print("MRFI Roots_den_eval: ",Roots_den_eval):

    # Generate monomials and coefficients using hash tables (already initialized at the start)
    # print("lambda_num: ",lambda_num):
    # print("lambda_den: ",lambda_den):
    
    for k from 1 to num_eqn do
        print("-----------------------------------"):
        temp := generate_monomials(Roots_num_eval[k], num_vars, Primes, vars):
        print("temp= ",temp):
        if temp = FAIL then return FAIL: end if:
        
        # Store in hash table with key k
        num_mono[k] := temp:
        print("MRFI num_mono[",k,"]:", num_mono[k]):
        print("MRFI lambda_num[",k,"]:", lambda_num[k]):
        
        coeff_num[k] := Zippel_Transpose_Vandermonde_solver(num_eval[k], terms_num[k], 
                                     Roots_num_eval[k], lambda_num[k], p):
        
    end do:

    # print("num_mono: ",num_mono):
    # print("coeff_num: ",coeff_num):
    # print("final_num: ",final_num):
    
    # Generate denominators using hash tables (already initialized at the start)
    if common_den_flag then
        temp := generate_monomials(Roots_den_eval[1], num_vars, Primes, vars):
        if temp = FAIL then return FAIL: end if:
        
        # Store common monomial list for first denominator
        den_mono[1] := temp:
        
        coeff_den[1] := Zippel_Transpose_Vandermonde_solver(den_eval[1], terms_den[1],
                                               Roots_den_eval[1], lambda_den[1], p):
        
        temp_den := construct_final_polynomial(coeff_den[1], den_mono[1]):
        
        # Set the same denominator for all equations
        for k from 1 to num_eqn do
            final_den[k] := temp_den:
            if k > 1 then
                den_mono[k] := den_mono[1]:
                coeff_den[k] := coeff_den[1]:
            end if:
        end do:
    else
        for k from 1 to num_eqn do
            temp := generate_monomials(Roots_den_eval[k], num_vars, Primes, vars):
            if temp = FAIL then return FAIL: end if:
            
            # Store in hash table with key k
            den_mono[k] := temp:
            
            coeff_den[k] := Zippel_Transpose_Vandermonde_solver(den_eval[k], terms_den[k],
                                               Roots_den_eval[k], lambda_den[k], p):
            

        end do:
    end if:

    for k from 1 to num_eqn do
        # print("k= ",k):
        # print("coeff_num[",k,"]",coeff_num[k]):
        # print("coeff_den[",k,"]",coeff_den[k]):
        # print("After normalizing"):

        # final_num[k] := construct_final_polynomial( coeff_num[k], num_mono[k] ):
        # final_den[k] := construct_final_polynomial( coeff_den[k], den_mono[k] ):
        # u:=1/lcoeff(final_den[k],vars) mod p:
        # final_num[k]:=u*final_num[k] mod p:
        # final_den[k]:=u*final_den[k] mod p:


        lcoeff( add( mon, mon in den_mono[k] ), vars, 'mon' );
        if not member( mon, den_mono[k], 'i' ) then error "bug in leading monomial" fi;
#         printf("LMON %d %a\n",i,mon);
        u := 1/coeff_den[k][-1] mod p;
        # u := 1/coeff_den[k][i] mod p;
        # printf(" LC(den)=%d  u=%d\n", coeff_den[k][-1], u);
        coeff_num[k]:=u*coeff_num[k] mod p:
        coeff_den[k]:=u*coeff_den[k] mod p:
# printf("CHECK  N := %a; DD := %a;\n", 
        iratrecon(add(coeff_num[k][i]*num_mono[k][i],i=1..nops(coeff_num[k])),p),
        iratrecon(add(coeff_den[k][i]*den_mono[k][i],i=1..nops(coeff_den[k])),p);
        print("coeff_num[",k,"]",coeff_num[k]):
        print("coeff_den[",k,"]",coeff_den[k]):
        final_num[k] := construct_final_polynomial( coeff_num[k], num_mono[k] ):
        final_den[k] := construct_final_polynomial( coeff_den[k], den_mono[k] ):
    end do:
    
    lprint("MRFI ========================================"):
    lprint("MRFI RECOVERY COMPLETE"):
    lprint("MRFI ========================================"):

    return final_num,final_den:

end proc: