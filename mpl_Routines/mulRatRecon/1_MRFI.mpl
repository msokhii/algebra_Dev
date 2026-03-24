MRFI := proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, Primes, direction, sigma_, num_eval, den_eval, u, mon,
          numerator_done, denominator_done, T, T_old,
          mqrfr_results, lin_sys, num_points_mqrfr,
          Numerators, Denominiators, deg_num, deg_den,
          lambda_num, lambda_den, terms_num, terms_den,
          R_num, R_den, Roots_num_eval, Roots_den_eval,
          num_mono, den_mono, coeff_num, coeff_den,
          final_num, final_den, r_, temp, common_den_flag,
          den_lc, den_lc_inv, bmea_done, r, temp_den,
          all_den_done, all_done, max_num_points_mqrfr,
          tempTable, ck, tempN, tempD, rTemp, alpha, Psi_alpha,
          Y, U, m, iter, iter2, mqrfr_list;

    lprint("MRFI  ========================================"):
    lprint("MRFI Starting MRFI"):
    lprint("MRFI Number of parameters:", num_vars):
    lprint("MRFI Number of equations:", num_eqn):
    lprint("MRFI ========================================"):

    r_ := rand(p):
    Primes := [seq(ithprime(i), i=1..num_vars)]:
    direction := [seq(r_(), i=1..num_vars-1)]:
    sigma_ := []:

    num_eval := Array(1..num_eqn):
    den_eval := Array(1..num_eqn):
    numerator_done := Array(1..num_eqn):
    denominator_done := Array(1..num_eqn):
    bmea_done := Array(1..num_eqn):
    num_points_mqrfr := Array(1..num_eqn):

    for i from 1 to num_eqn do
        num_eval[i] := []:
        den_eval[i] := []:
        numerator_done[i] := false:
        denominator_done[i] := false:
        bmea_done[i] := false:
        num_points_mqrfr[i] := 0:
    end do:

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

    lprint("MRFI numerator_done: ", [seq(numerator_done[i], i=1..num_eqn)]):
    lprint("MRFI denominator_done: ", [seq(denominator_done[i], i=1..num_eqn)]):

    T := 4:
    T_old := 1:

    mqrfr_results, lin_sys := NDSA(B, [seq(1, i=1..num_vars)], direction, num_vars, p, T, num_eqn):
    print(mqrfr_results):
    print(lin_sys):

    Numerators := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    Denominiators := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:

    for k from 1 to num_eqn do
        num_eval[k] := [eval(Numerators[k], x=1) mod p]:
        den_eval[k] := [eval(Denominiators[k], x=1) mod p]:
    end do:

    deg_num := [seq(degree(Numerators[i], x), i=1..nops(Numerators))]:
    deg_den := [seq(degree(Denominiators[i], x), i=1..nops(Denominiators))]:

    print(deg_num):
    print(deg_den):
    print(T_old):
    nops(deg_num);

    for i from 1 to nops(deg_den) do
        num_points_mqrfr[i] := deg_num[i] + deg_den[i] + 2:
    end do:

    max_num_points_mqrfr := max(op(deg_num)) + max(op(deg_den)) + 2:
    print("MRFI num_points_mqrfr: ", [seq(num_points_mqrfr[i], i=1..num_eqn)]):
    print("MRFI max_num_points_mqrfr: ", max_num_points_mqrfr):

    common_den_flag := true:
    for k from 2 to num_eqn do
        if Denominiators[k] <> Denominiators[1] then
            common_den_flag := false:
            break:
        end if:
    end do:

    print("______________________________________________________________________________"):
    print("in main evaluation loop of MRFI"):

    while true do
        lprint("MRFI T_old=", T_old):
        lprint("MRFI T=", T):

        tempTable := table():
        ck := 1:
        
        for j from T_old to 2*T-1 do
            sigma_ := [op(sigma_), [seq(Primes[i]^j mod p, i=1..nops(Primes))]]:
            print(sigma_):

            if j <= nops(deg_num) then
                tempN := deg_num[j]:
            else
                tempN := max(op(deg_num)):
            end if:
            print(tempN):

            if j <= nops(deg_den) then
                tempD := deg_den[j]:
            else
                tempD := max(op(deg_den)):
            end if:
            print(tempD):

            rTemp := rand(p):
            alpha := [seq(rTemp(), iter=1..tempN+tempD+1)]:
            print(alpha):

            Psi_alpha := get_point_on_affine_line_2(num_vars, alpha, direction, sigma_[j], p,nops(alpha)):
            print(Psi_alpha):
            Y := [seq(B(Psi_alpha[i], p), i=1..nops(Psi_alpha))]:
            U := Interp(alpha, Y, x) mod p:
            m := Expand(product(x-alpha[iter2], iter2=1..nops(alpha))) mod p:

            tempTable[ck] := RatRecon(U, m, x, tempN, tempD) mod p:
            ck := ck + 1:
        end do:
        tempTable := convert(tempTable,list):
        print(tempTable):
        print(nops(tempTable)):

        quit;
        mqrfr_list := [seq(tempTable[i], i=1..ck-1)]:
        Numerators := [seq(mqrfr_list[i][1], i=1..nops(mqrfr_list))]:
        Denominiators := [seq(mqrfr_list[i][2], i=1..nops(mqrfr_list))]:

        print("______________________________________________________________________________"):

        for k from 1 to num_eqn do
            if k <= nops(Numerators) and not numerator_done[k] then
                num_eval[k] := [op(num_eval[k]), eval(Numerators[k], x=sigma_[j][1]) mod p]:
            end if:
            if k <= nops(Denominiators) and not denominator_done[k] then
                den_eval[k] := [op(den_eval[k]), eval(Denominiators[k], x=sigma_[j][1]) mod p]:
            end if:
        end do:

        print("MRFI num_eval: ", [seq(num_eval[i], i=1..num_eqn)]):
        print("MRFI den_eval: ", [seq(den_eval[i], i=1..num_eqn)]):
        print("______________________________________________________________________________"):
        print("numerator_done: ", [seq(numerator_done[i], i=1..num_eqn)]):
        print("denominator_done: ", [seq(denominator_done[i], i=1..num_eqn)]):
        print("--------------------------------------------------------------------------------"):

        all_done := true:

        for k from 1 to num_eqn do
            lprint("MRFI k=", k):
            lprint(cat("MRFI numerator_done[", k, "]="), numerator_done[k]):

            if numerator_done[k] then
                print("MRFI Skipping BMEA for numerator component ", k, " as already done"):
                next:
            end if:

            lambda_num[k] := BMEA(num_eval[k], p, Z):
            terms_num[k] := degree(lambda_num[k], Z):
            R_num[k] := Roots(lambda_num[k]) mod p:

            lprint("MRFI lambda_num: ", lambda_num[k]):
            lprint("MRFI terms_num: ", terms_num[k]):
            lprint("MRFI R_num: ", R_num[k]):
            print("--------------------------------------------------------------------------------"):
            print("Checking termination condition for numerator"):
            print(cat("nops(R_num[", k, "]): "), nops(R_num[k])):
            print(cat("R_num[", k, "]: "), R_num[k]):
            print(cat("terms_num[", k, "]: "), terms_num[k]):
            print("T: ", T):

            if R_num[k] = [] then
                print("MRFI: Empty roots list for numerator component ", k):
                numerator_done[k] := false:
                all_done := false:
                next:
            end if:

            if nops(R_num[k]) > 0 and R_num[k][1][1] = 0 then
                R_num[k] := remove(x -> x = [0,1], R_num[k]):
                terms_num[k] := terms_num[k] - 1:
            end if:

            if nops(R_num[k]) = terms_num[k] and terms_num[k] <= T then
                print("MRFI Numerator component ", k, " recovered!"):
                numerator_done[k] := true:
            else
                all_done := false:
            end if:

            print("____________________________________________________________________________"):
        end do:

        print("===================================================================================="):
        print("MRFI processing denominators now"):
        lprint("MRFI den_eval: ", [seq(den_eval[i], i=1..num_eqn)]):

        all_den_done := true:
        lprint("MRFI common_den_flag: ", common_den_flag):

        if common_den_flag then
            lprint("In common_den_flag"):

            if not denominator_done[1] then
                lambda_den[1] := BMEA(den_eval[1], p, Z):
                terms_den[1] := degree(lambda_den[1], Z):
                R_den[1] := Roots(lambda_den[1]) mod p:

                lprint("MRFI lambda_den: ", lambda_den[1]):
                lprint("MRFI terms_den: ", terms_den[1]):
                lprint("MRFI R_den: ", R_den[1]):

                if nops(R_den[1]) > 0 and R_den[1][1][1] = 0 then
                    R_den[1] := remove(x -> x = [0,1], R_den[1]):
                    terms_den[1] := terms_den[1] - 1:
                end if:

                if nops(R_den[1]) = terms_den[1] and terms_den[1] < T then
                    for k from 1 to num_eqn do
                        denominator_done[k] := true:
                    end do:
                else
                    all_den_done := false:
                end if:
            end if:
        else
            for k from 1 to num_eqn do
                if denominator_done[k] then
                    print("MRFI Skipping BMEA for denominator component ", k, " as already done"):
                    next:
                end if:

                lambda_den[k] := BMEA(den_eval[k], p, Z):
                terms_den[k] := degree(lambda_den[k], Z):
                R_den[k] := Roots(lambda_den[k]) mod p:

                print("--------------------------------------------------------------------------------"):
                print("Checking termination condition for denominator"):
                print(cat("terms_den[", k, "]: "), terms_den[k]):
                print("T: ", T):

                if nops(R_den[k]) > 0 and R_den[k][1][1] = 0 then
                    R_den[k] := remove(x -> x = [0,1], R_den[k]):
                    terms_den[k] := terms_den[k] - 1:
                end if:

                if nops(R_den[k]) = terms_den[k] and terms_den[k] < T then
                    print("MRFI Denominator component ", k, " recovered!"):
                    denominator_done[k] := true:
                else
                    all_den_done := false:
                end if:
            end do:
        end if:

        print("--------------------------------------------------------------------------------"):
        for i from 1 to num_eqn do
            bmea_done[i] := numerator_done[i] and denominator_done[i]:
        end do:

        print("BMEA done status: ", [seq(bmea_done[i], i=1..num_eqn)]):

        all_done := true:
        for i from 1 to num_eqn do
            all_done := all_done and bmea_done[i]:
        end do:

        print("All done status: ", all_done):

        if all_done then
            lprint("MRFI All components recovered!"):
            break:
        end if:

        print("numerator_done: ", [seq(numerator_done[i], i=1..num_eqn)]):
        print("denominator_done: ", [seq(denominator_done[i], i=1..num_eqn)]):
        print("R_num: ", R_num):
        print("R_den: ", R_den):
        print("terms_num: ", terms_num):
        print("terms_den: ", terms_den):
        print("lambda_den: ", lambda_den):

        T_old := 2*T:
        T := T*2:
        print("______________________________________________________________________________"):

        if T > 2^6 then
            break:
        end if:
    end do:

    Roots_num_eval := [seq([seq(r[1], r in R_num[k])], k=1..num_eqn)]:
    print("MRFI Roots_num_eval: ", Roots_num_eval):

    if common_den_flag then
        Roots_den_eval := [[seq(r[1], r in R_den[1])]]:
        for k from 2 to num_eqn do
            Roots_den_eval := [op(Roots_den_eval), Roots_den_eval[1]]:
        end do:
    else
        Roots_den_eval := [seq([seq(r[1], r in R_den[k])], k=1..num_eqn)]:
    end if:

    print("MRFI Roots_den_eval: ", Roots_den_eval):

    for k from 1 to num_eqn do
        print("-----------------------------------"):
        temp := generate_monomials(Roots_num_eval[k], num_vars, Primes, vars):
        print("temp= ", temp):
        if temp = FAIL then
            return FAIL:
        end if:

        num_mono[k] := temp:
        print(cat("MRFI num_mono[", k, "]:"), num_mono[k]):
        print(cat("MRFI lambda_num[", k, "]:"), lambda_num[k]):

        coeff_num[k] := Zippel_Transpose_Vandermonde_solver(
            num_eval[k], terms_num[k], Roots_num_eval[k], lambda_num[k], p
        ):
    end do:

    if common_den_flag then
        temp := generate_monomials(Roots_den_eval[1], num_vars, Primes, vars):
        if temp = FAIL then
            return FAIL:
        end if:

        den_mono[1] := temp:
        coeff_den[1] := Zippel_Transpose_Vandermonde_solver(
            den_eval[1], terms_den[1], Roots_den_eval[1], lambda_den[1], p
        ):
        temp_den := construct_final_polynomial(coeff_den[1], den_mono[1]):

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
            if temp = FAIL then
                return FAIL:
            end if:

            den_mono[k] := temp:
            coeff_den[k] := Zippel_Transpose_Vandermonde_solver(
                den_eval[k], terms_den[k], Roots_den_eval[k], lambda_den[k], p
            ):
        end do:
    end if:

    for k from 1 to num_eqn do
        lcoeff(add(mon, mon in den_mono[k]), vars, 'mon'):
        if not member(mon, den_mono[k], 'i') then
            error "bug in leading monomial":
        end if:

        u := 1/coeff_den[k][-1] mod p:
        coeff_num[k] := u*coeff_num[k] mod p:
        coeff_den[k] := u*coeff_den[k] mod p:

        printf(
            "CHECK  N := %a; DD := %a;\n",
            iratrecon(add(coeff_num[k][i]*num_mono[k][i], i=1..nops(coeff_num[k])), p),
            iratrecon(add(coeff_den[k][i]*den_mono[k][i], i=1..nops(coeff_den[k])), p)
        ):

        print(cat("coeff_num[", k, "]"), coeff_num[k]):
        print(cat("coeff_den[", k, "]"), coeff_den[k]):

        final_num[k] := construct_final_polynomial(coeff_num[k], num_mono[k]):
        final_den[k] := construct_final_polynomial(coeff_den[k], den_mono[k]):
    end do:

    lprint("MRFI ========================================"):
    lprint("MRFI RECOVERY COMPLETE"):
    lprint("MRFI ========================================"):

    return final_num, final_den:
end proc:
