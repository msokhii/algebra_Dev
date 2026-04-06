rrMRFI:= proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, s, u, temp,
          Primes, direction, sigma_j,
          numerator_done, denominator_done, bmea_done, all_done,
          Tinit, Tcur, jDone,
          mqrfr_results, lin_sys,
          Numerators, Denominators, deg_num, deg_den,
          sampleCounts, m_i, mMax,
          lambda_num, lambda_den, terms_num, terms_den,
          R_num, R_den,
          coeff_num, coeff_den,
          common_den_flag,
          tempG, ratReconVal,
          num_eval, den_eval,
          alphaVal, Psi_alpha, BBvals, Y, interpVal, rr,
          modulusTable, modulusPoly,
          r_, initialPoint,
          Roots_num_eval, Roots_den_eval,
          num_mono, den_mono, final_num, final_den,
          tmpNum, tmpDen,
          maxDoublings, doublingCount;

    lprint("RR MRFI ========================================"):
    lprint("RR MRFI Starting"):
    lprint("num_vars =", num_vars):
    lprint("num_eqn  =", num_eqn):
    lprint("=============================================================="):

    r_ := rand(p):
    tempG := rand(p):

    Primes := [seq(ithprime(i), i=1..num_vars)]:
    direction := [seq(r_(), i=1..num_vars-1)]:
    initialPoint := [seq(1, i=1..num_vars)]:

    numerator_done   := [seq(false, i=1..num_eqn)]:
    denominator_done := [seq(false, i=1..num_eqn)]:
    bmea_done        := [seq(false, i=1..num_eqn)]:

    lambda_num := table():
    terms_num  := table():
    R_num      := table():
    lambda_den := table():
    terms_den  := table():
    R_den      := table():
    coeff_num  := table():
    coeff_den  := table():
    ratReconVal := table():
    num_mono   := table():
    den_mono   := table():
    final_num  := table():
    final_den  := table():

    for i from 1 to num_eqn do
        lambda_num[i] := []:
        terms_num[i]  := 0:
        R_num[i]      := []:
        lambda_den[i] := []:
        terms_den[i]  := 0:
        R_den[i]      := []:
        coeff_num[i]  := 0:
        coeff_den[i]  := 0:
        ratReconVal[i] := table():
        num_mono[i]   := []:
        den_mono[i]   := []:
        final_num[i]  := 0:
        final_den[i]  := 0:
    end do:

    # ============================================================
    # Phase 1: initial degree discovery from NDSA
    # ============================================================
    Tinit := 4:
    mqrfr_results, Tcur, lin_sys := NDSA(B, initialPoint, direction, num_vars, p, Tinit, num_eqn):

    Numerators   := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    Denominators := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:

    deg_num := [seq(degree(Numerators[i], x), i=1..nops(Numerators))]:
    deg_den := [seq(degree(Denominators[i], x), i=1..nops(Denominators))]:

    sampleCounts := [seq(deg_num[i] + deg_den[i] + 1, i=1..num_eqn)]:
    mMax := max(op(sampleCounts)):

    print("deg_num      =", deg_num):
    print("deg_den      =", deg_den):
    print("sampleCounts =", sampleCounts):
    print("mMax         =", mMax):
    print("T from NDSA  =", Tcur):

    common_den_flag := true:
    for k from 2 to num_eqn do
        if Denominators[k] <> Denominators[1] then
            common_den_flag := false:
            break:
        end if:
    end do:
    print("common_den_flag =", common_den_flag):

    # ============================================================
    # Seed num_eval / den_eval with x = 1 values from NDSA output
    # ============================================================
    num_eval := [seq([eval(Numerators[k], x=1) mod p], k=1..num_eqn)]:
    den_eval := [seq([eval(Denominators[k], x=1) mod p], k=1..num_eqn)]:

    # ============================================================
    # Phase 2: adaptive outer loop on T
    # ============================================================
    jDone := 0:
    all_done := false:
    maxDoublings := 20:
    doublingCount := 0:

    while not all_done do
        print("=============================================================="):
        print("Current Tcur =", Tcur):
        print("Current jDone =", jDone):
        print("Will process j from", jDone+1, "to", 2*Tcur-1):

        # --------------------------------------------------------
        # Only process NEW j-values. Do not recompute old ones.
        # --------------------------------------------------------
        for j from jDone+1 to 2*Tcur-1 do
            sigma_j := [seq(Primes[k]^j mod p, k=1..nops(Primes))]:

            # One shared alpha list of length mMax
            alphaVal := [seq(tempG(), s=1..mMax)]:

            # Build prefix modulus polynomials:
            # modulusTable[s] = product_{t=1..s} (x - alphaVal[t])
            modulusTable := table():
            modulusPoly := 1:
            for s from 1 to mMax do
                modulusPoly := Expand(modulusPoly*(x-alphaVal[s])) mod p:
                modulusTable[s] := modulusPoly:
            end do:

            # One shared affine-line sample set
            Psi_alpha := get_point_on_affine_line(num_vars, alphaVal, direction, sigma_j, p, mMax):

            # One shared BB batch
            BBvals := [seq(B(Psi_alpha[s], p), s=1..mMax)]:

            # Reconstruct each component using only its needed prefix length
            for i from 1 to num_eqn do
                m_i := sampleCounts[i]:

                Y := [seq(BBvals[s][i], s=1..m_i)]:

                interpVal := cppNewtonInterp([op(1..m_i, alphaVal)], Y, x, p):
                rr := cppRR(interpVal, modulusTable[m_i], x, deg_num[i], deg_den[i],p):

                ratReconVal[i][j] := rr:

                tmpNum := Eval(numer(rr), x=sigma_j[1]) mod p:
                tmpDen := Eval(denom(rr), x=sigma_j[1]) mod p:

                num_eval[i] := [op(num_eval[i]), tmpNum]:
                den_eval[i] := [op(den_eval[i]), tmpDen]:
            end do:
        end do:

        jDone := 2*Tcur - 1:

        # --------------------------------------------------------
        # Reset done flags before BMEA check on accumulated data
        # --------------------------------------------------------
        for k from 1 to num_eqn do
            numerator_done[k]   := false:
            denominator_done[k] := false:
            bmea_done[k]        := false:
        end do:

        print("MRFI2_global num_eval:", num_eval):
        print("MRFI2_global den_eval:", den_eval):
        print("______________________________________________________________________________"):
        print("numerator_done: ", numerator_done):
        print("denominator_done: ", denominator_done):
        print("--------------------------------------------------------------------------------"):

        all_done := true:

        # ============================================================
        # Process numerators
        # ============================================================
        for k from 1 to num_eqn do
            lprint("MRFI k=", k):
            lprint("MRFI numerator_done[", k, "]=", numerator_done[k]):

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
            print("nops(R_num[", k, "]): ", nops(R_num[k])):
            print("R_num[", k, "]: ", R_num[k]):
            print("terms_num[", k, "]: ", terms_num[k]):
            print("Tcur: ", Tcur):

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

            if nops(R_num[k]) = terms_num[k] and terms_num[k] <= Tcur then
                print("MRFI Numerator component ", k, " recovered!"):
                numerator_done[k] := true:
            else
                all_done := false:
            end if:

            print("____________________________________________________________________________"):
        end do:

        print("===================================================================================="):
        print("MRFI processing denominators now"):

        lprint("MRFI den_eval: ", den_eval):
        lprint("MRFI common_den_flag: ", common_den_flag):

        # ============================================================
        # Process denominators
        # ============================================================
        if common_den_flag then
            lprint("In common_den_flag"):

            if denominator_done[1] then
                print("MRFI Skipping common denominator BMEA as already done"):
            else
                lambda_den[1] := BMEA(den_eval[1], p, Z):
                terms_den[1] := degree(lambda_den[1], Z):
                R_den[1] := Roots(lambda_den[1]) mod p:

                lprint("MRFI lambda_den: ", lambda_den[1]):
                lprint("MRFI terms_den: ", terms_den[1]):
                lprint("MRFI R_den: ", R_den[1]):

                if R_den[1] = [] then
                    print("MRFI: Empty roots list for common denominator"):
                    all_done := false:
                else
                    if nops(R_den[1]) > 0 and R_den[1][1][1] = 0 then
                        R_den[1] := remove(x -> x = [0,1], R_den[1]):
                        terms_den[1] := terms_den[1] - 1:
                    end if:

                    if nops(R_den[1]) = terms_den[1] and terms_den[1] < Tcur then
                        for k from 1 to num_eqn do
                            denominator_done[k] := true:
                        end do:
                    else
                        all_done := false:
                    end if:
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

                lprint("MRFI lambda_den[", k, "]: ", lambda_den[k]):
                lprint("MRFI terms_den[", k, "]: ", terms_den[k]):
                lprint("MRFI R_den[", k, "]: ", R_den[k]):

                print("--------------------------------------------------------------------------------"):
                print("Checking termination condition for denominator"):
                print("terms_den[", k, "]: ", terms_den[k]):
                print("Tcur: ", Tcur):

                if R_den[k] = [] then
                    print("MRFI: Empty roots list for denominator component ", k):
                    denominator_done[k] := false:
                    all_done := false:
                    next:
                end if:

                if nops(R_den[k]) > 0 and R_den[k][1][1] = 0 then
                    R_den[k] := remove(x -> x = [0,1], R_den[k]):
                    terms_den[k] := terms_den[k] - 1:
                end if:

                if nops(R_den[k]) = terms_den[k] and terms_den[k] < Tcur then
                    print("MRFI Denominator component ", k, " recovered!"):
                    denominator_done[k] := true:
                else
                    all_done := false:
                end if:
            end do:
        end if:

        print(denominator_done):
        print(numerator_done):
        print("--------------------------------------------------------------------------------"):

        for i from 1 to num_eqn do
            bmea_done[i] := numerator_done[i] and denominator_done[i]:
        end do:

        print("BMEA done status: ", bmea_done):

        all_done := true:
        for i from 1 to num_eqn do
            all_done := all_done and bmea_done[i]:
        end do:

        print("All done status: ", all_done):
        if all_done then
            lprint("MRFI All components recovered!"):
        end if:

        print("numerator_done: ", numerator_done):
        print("denominator_done: ", denominator_done):
        print("R_num: ", R_num):
        print("R_den: ", R_den):
        print("terms_num: ", terms_num):
        print("terms_den: ", terms_den):
        print("lambda_num: ", lambda_num):
        print("lambda_den: ", lambda_den):

        if not all_done then
            doublingCount := doublingCount + 1:
            if doublingCount > maxDoublings then
                error "Exceeded maximum number of T doublings without completion":
            end if:

            Tcur := 2*Tcur:
            print("Doubling Tcur to", Tcur):
        end if:
    end do:

    # ============================================================
    # Convert root structures to evaluation root lists
    # ============================================================
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

    # ============================================================
    # Generate numerator monomials and coefficients
    # ============================================================
    for k from 1 to num_eqn do
        print("-----------------------------------"):
        temp := generate_monomials(Roots_num_eval[k], num_vars, Primes, vars):
        print("temp= ", temp):
        if temp = FAIL then
            return FAIL:
        end if:

        num_mono[k] := temp:
        print("MRFI num_mono[", k, "]:", num_mono[k]):
        print("MRFI lambda_num[", k, "]:", lambda_num[k]):

        coeff_num[k] := Zippel_Transpose_Vandermonde_solver(
                            num_eval[k],
                            terms_num[k],
                            Roots_num_eval[k],
                            lambda_num[k],
                            p
                        ):
    end do:

    # ============================================================
    # Generate denominator monomials and coefficients
    # ============================================================
    if common_den_flag then
        temp := generate_monomials(Roots_den_eval[1], num_vars, Primes, vars):
        if temp = FAIL then
            return FAIL:
        end if:

        den_mono[1] := temp:

        coeff_den[1] := Zippel_Transpose_Vandermonde_solver(
                            den_eval[1],
                            terms_den[1],
                            Roots_den_eval[1],
                            lambda_den[1],
                            p
                        ):

        temp := construct_final_polynomial(coeff_den[1], den_mono[1]):

        for k from 1 to num_eqn do
            final_den[k] := temp:
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
                                den_eval[k],
                                terms_den[k],
                                Roots_den_eval[k],
                                lambda_den[k],
                                p
                            ):
        end do:
    end if:

    # ============================================================
    # Normalize and construct final numerator / denominator
    # ============================================================
    for k from 1 to num_eqn do
        lcoeff(add(mon, mon in den_mono[k]), vars, 'mon'):
        if not member(mon, den_mono[k], 'i') then
            error "bug in leading monomial":
        end if:

        u := 1/coeff_den[k][-1] mod p:

        coeff_num[k] := u*coeff_num[k] mod p:
        coeff_den[k] := u*coeff_den[k] mod p:

        print("coeff_num[", k, "]", coeff_num[k]):
        print("coeff_den[", k, "]", coeff_den[k]):

        final_num[k] := construct_final_polynomial(coeff_num[k], num_mono[k]):
        final_den[k] := construct_final_polynomial(coeff_den[k], den_mono[k]):
    end do:

    lprint("MRFI ========================================"):
    lprint("MRFI RECOVERY COMPLETE"):
    lprint("MRFI ========================================"):
    print("Final Tcur =", Tcur):
    print("Final jDone =", jDone):
    print("Main-loop probe estimate =", jDone * mMax):

    return final_num, final_den:
end proc: