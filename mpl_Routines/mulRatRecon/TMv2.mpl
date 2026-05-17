kernelopts(numcpus=1):
read "./mapleWrapperv2.mpl":

with(LinearAlgebra):
with(IntegerRelations):

get_eqn := proc(Sys, vars)
    option remember;
    return solve(Sys, vars);
end proc:

reording := proc(unordered_soln, num_eqn)
    local component1, ordered_soln, i;
    ordered_soln := [seq(0, i=1..num_eqn)];
    for i from 1 to num_eqn do
        component1 := get_component(op(1, unordered_soln[i]));
        ordered_soln[component1] := unordered_soln[i];
    end do;
    return ordered_soln;
end proc:

get_component := proc(expression)
    local temp_var;
    temp_var := convert(expression, string);
    return parse(temp_var[2..length(temp_var)]);
end proc:

get_u := proc(M, col, alpha, p)
    local F, U, i;
    F := [seq(convert(M[..,i], list), i=1..col)];
    U := [seq(cppNewtonInterp(alpha,F[i],x,p),i=1..col)];
    return U;
end proc:

Deterministic_get_u := proc(M, col, alpha, p)
    local F, U;
    F := convert(M[..,col][..numelems(alpha)], list);
    U := cppNewtonInterp(alpha,F,x,p);
    return U;
end proc:

FFGE := proc(A::Matrix, b::Vector, Y::list(name))
local n, m, B, mu, i, j, k, det, x, y, num, r, numterms, f, g, h, mons,
      elim_num_max, elim_post_max, backsub_N_max, y_terms, f_terms, g_terms;
    global ffge_stats;
    numterms := proc(f) if f=0 then 0 elif type(f,`+`) then nops(f) else 1 fi end;
    n,m := op(1,A);
    if n<>m then error "Matrix must be square" fi;
    m := op(1,b);
    if m<>n then error "Matrix and vector must have the same dimension" fi;
    B := <A|b>;
    mu := 1:
    det := 1;
    elim_num_max  := 0:
    elim_post_max := 0:
    for k to n-1 do
        i := k;
        while i<=n and B[i,k]=0 do i := i+1; od;
        if i>n then return 0 fi;
        if i>k then 
            for j from k to n+1 do B[i,j],B[k,j] := B[k,j],B[i,j] od;
            det := -det;
        fi;
        for i from k+1 to n do
            for j from k+1 to n+1 do
                num := expand(B[k,k]*B[i,j]-B[i,k]*B[k,j]);
                divide(num,mu,evaln(B[i,j]));
                if i=k+1 and j=k+1 then
                    #lprint(i,numterms(num),numterms(B[i,j]));
                    if numterms(num)      > elim_num_max  then elim_num_max  := numterms(num)      end if;
                    if numterms(B[i,j])   > elim_post_max then elim_post_max := numterms(B[i,j])   end if;
                fi;
            od;
            B[i,k] := 0;
        od;
        mu := B[k,k];
    od;
    det := det*B[n,n];
    #printf("#det=%d\n",numterms(det));
    # Lipson's back substitution
    y := Vector(n);
    y[n] := B[n,n+1];
    #printf("#y[%d]=%d\n",n,numterms(B[n,n+1]));
    backsub_N_max := 0:
    for i from n-1 by -1 to 1 do
        num := expand(B[i,n+1]*B[n,n]-add(B[i,j]*y[j],j=i+1..n));
        divide(num,B[i,i],evaln(y[i]));
        #printf("#N[%d]=%d  #y[%d]=%d\n",i,numterms(num),i,numterms(y[i]));
        if numterms(num) > backsub_N_max then backsub_N_max := numterms(num) end if;
    od;
    f := Vector(n);
    g := Vector(n);
    for i from 1 to n do
        h := gcd(y[i],B[n,n],evaln(f[i]),evaln(g[i]));
        #printf("#f[%d]=%d #g[%d]=%d #h=%d",i,numterms(f[i]),
             #i,numterms(g[i]),numterms(h));
        #printf("   deg(f[%d])=%d deg(g[%d])=%d\n",i,degree(f[i]),i,degree(g[i]));
    od;
    r := {seq(Y[i]=ithprime(i),i=1..nops(Y))};
    for i from 1 to n do
        coeffs(f[i],Y,'mons');
        mons := subs(r,[mons]);
        m := max(op(mons));
        #printf("f[%d]: max m[%d] = %d = %a\n",i,i,m,ifactor(m));
        coeffs(g[i],Y,'mons');
        mons := subs(r,[mons]);
        m := max(op(mons));
        #printf("g[%d]: max m[%d] = %d = %a\n",i,i,m,ifactor(m));
    od;

    # Publish stats for the driver.
    y_terms := [seq(numterms(y[i]), i=1..n)]:
    f_terms := [seq(numterms(f[i]), i=1..n)]:
    g_terms := [seq(numterms(g[i]), i=1..n)]:
    ffge_stats := table():
    ffge_stats["elim_num_max"]  := elim_num_max:
    ffge_stats["elim_post_max"] := elim_post_max:
    ffge_stats["backsub_N_max"] := backsub_N_max:
    ffge_stats["y_terms"]       := y_terms:
    ffge_stats["det_terms"]     := numterms(det):
    ffge_stats["f_terms"]       := f_terms:
    ffge_stats["g_terms"]       := g_terms:

    return det,f,g;
end:

Constuct_Sys_Blackbox := proc(Sys, Vars, params)
    local Lin_BB, L;
    L := GenerateMatrix(Sys, Vars, augmented=true):
    Lin_BB := proc(point_::list(integer), p::prime)
        local A, T, subs_values, num_eqn, soln, t0, t_avg, i;
        uses LinearAlgebra:-Modular;
        global counter, bb_calls_ndsa, bb_calls_mrfi,
               t_ndsa_total, t_mrfi_total, bb_phase, bb_bench_calls;
        counter := counter + 1;
        if bb_phase = "NDSA" then
            bb_calls_ndsa := bb_calls_ndsa + 1;
        else
            bb_calls_mrfi := bb_calls_mrfi + 1;
        end if;
        subs_values := zip((par, pnt) -> par = pnt, params, point_):
        num_eqn := numelems(Vars):

        # Time the BB over bb_bench_calls iterations (last iter's A is the result).
        t0 := time();
        for i from 1 to bb_bench_calls do
            A := Mod(p, L, subs_values, integer[8]):
            T := traperror(LinearSolve(p, A, 1)):
        end do;
        t_avg := evalf((time() - t0) / bb_bench_calls);

        if bb_phase = "NDSA" then
            t_ndsa_total := t_ndsa_total + t_avg;
        else
            t_mrfi_total := t_mrfi_total + t_avg;
        end if;

        if T = "matrix is singular" then return FAIL fi:
        soln := convert(A[1..num_eqn, num_eqn+1], list):
        return soln;
    end:
end:

get_point_on_affine_line := proc(num_var::posint, alpha::list, beta_::list,
                                 sigma_::list, p::prime, T::posint)
    option remember;
    local psi, nv, np, i;
    global num_lines;
    num_lines := num_lines + 1;
    if num_var < 2 then error "num_var must be >= 2, got %1", num_var; end if;
    if nops(alpha) < T then error "Insufficient alpha values: need %1, got %2",
                                  T, nops(alpha); end if;
    ASSERT(nops(beta_)  = num_var - 1, "beta_ must have num_var-1 elements");
    ASSERT(nops(sigma_) >= num_var,    "sigma_ must have >= num_var elements");
    for np from 1 to T do
        psi[np][1] := alpha[np];
        for nv from 2 to num_var do
            psi[np][nv] := beta_[nv-1]*alpha[np]
                         - beta_[nv-1]*sigma_[1]
                         + sigma_[nv] mod p;
        end do;
    end do;
    return [seq(convert(psi[i], list), i=1..T)];
end proc:

MQRFR := proc(r0::polynom, r1::polynom, t0::integer, t1::integer, p::prime)
    local r, t, q, i, f, g, qmax, lcg;
    r[0] := r0; r[1] := r1; t[0] := t0; t[1] := t1;
    f := r0; g := t1; qmax := 0; i := 1;
    while r[i] <> 0 do
        q[i] := Quo(r[i-1], r[i], x, 'r[i+1]') mod p;
        if degree(q[i], x) > qmax then
            qmax := degree(q[i], x);
            f := r[i]; g := t[i];
        end if;
        t[i+1] := Expand(t[i-1] - q[i]*t[i]) mod p;
        i := i + 1;
        if qmax <= 1 or gcd(f,g) <> 1 or g = 0 then FAIL; end if;
    end do;
    lcg := lcoeff(g);
    return f/lcg mod p, g/lcg mod p, qmax, lcg;
end proc:

BMEA := proc(v::list, p::posint, Z::name)
    # Wraps the C++ Berlekamp-Massey via cppBMM, lifting the returned
    # coefficient list [Lambda_0, Lambda_1, ..., Lambda_deg] into a polynomial
    # in Z. Empty / all-zero input returns 1 so Roots(...) and degree(...)
    # downstream stay well-behaved.
    local L, d, i;
    if numelems(v) = 0 then return 1; end if;
    L := cppBMM(v, p);
    if L = [] then return 1; end if;
    d := nops(L) - 1;
    return add(L[i+1]*Z^i, i=0..d);
end:

generate_monomials := proc(roots_, num_var, prime_points, vars)
    local m, mm, i, j, counter_, M_, rem;
    M_ := Vector(numelems(roots_), 0);
    for i from 1 to numelems(roots_) do
        if roots_[i] = 0 then return FAIL; end if;
        mm := roots_[i];
        m := 1;
        for j from 1 to numelems(prime_points) do
            counter_ := 0;
            while mm mod prime_points[j] = 0 do
                mm := iquo(mm, prime_points[j], 'rem');
                counter_ := counter_ + 1;
            end do;
            m := m*vars[j]^counter_;
        end do;
        M_[i] := m;
        if mm <> 1 then
            print("Warning: residue mm=", mm, " (should be 1) for root[", i, "]");
            return FAIL;
        end if;
    end do;
    return convert(M_, list);
end proc:

Zippel_Transpose_Vandermonde_solver := proc(y::list, terms::integer,
                                            roots_::list, lambda_::polynom,
                                            p::integer)
    # Wraps the C++ transposed Vandermonde solve via cppVS. lambda_ is unused:
    # cppVS builds the master polynomial M(x) = prod(x-roots_[i]) internally.
    # shift=1 folds in the trailing 1/roots_[i] factor (the three MRFI call
    # sites no longer apply that division as a post-step).
    if terms = 0 then return []; end if;
    return cppVS(y[1..terms], roots_[1..terms], p, 1);
end proc:

construct_final_polynomial := proc(coeff_, Monomials)
    local i, f;
    f := 0;
    for i from 1 to numelems(coeff_) do
        f := f + coeff_[i]*Monomials[i];
    end do;
    return f;
end proc:

RVR := proc(L, p, field)
    local M, E, H, B, i, rvr, num_coeff, a, q, g, ok, qX;
    M := p; E := ceil(sqrt(p)); num_coeff := numelems(L);
    H := L; H := [op(H), E];
    B := Matrix(num_coeff, num_coeff, shape=identity);
    B := p*B;
    B := <B|Vector([seq(0, i=1..num_coeff)])>;
    B := <B; H>;
    rvr := LLL(B);
    if field = "Z" then return rvr[1, 1..-2];
    elif field = "Q" then
        qX := rvr[1, num_coeff+1];
        a := Vector(num_coeff, i -> rvr[1, i]);
        q := qX/E;
        g := igcd(q, seq(a[i], i=1..num_coeff));
        q := iquo(q, g);
        a := Vector(num_coeff, i -> iquo(a[i], g));
        if q < 0 then
            q := -q;
            a := Vector(num_coeff, i -> -a[i]);
        end if;
        ok := true;
        for i to num_coeff do
            if irem(a[i] - q*L[i], M) <> 0 then ok := false; break; end if;
        end do;
        if not ok then error "verification failed"; end if;
        return [seq(a[i]/q, i=1..num_coeff)];
    end if;
end proc:

NDSA := proc(B, sigma_, beta_, num_var, p, num_points, num_eqn)
    local correct_degree, T, alpha, m, Psi_alpha, Y, u, dq, i, r,
          lin_sys, temp, result, count, M, row, col, DQ, MQRFR_done,
          t_helper;
    global t_cppNewton_total, t_mqrfr_total;
    print("In NDSA");
    MQRFR_done := [seq(false, i=1..num_eqn)]:
    correct_degree := false:
    lin_sys := false:
    T := num_points;
    temp := []:
    result := [seq([], i=1..num_eqn)]:
    count := 0:
    while not correct_degree do
        count := count + 1;
        print("NDSA: T = ", T);
        r := rand(p);
        alpha := [seq(i mod p, i=1..T)]:
        m := Expand(product(x - alpha[j], j=1..T)) mod p:
        Psi_alpha := get_point_on_affine_line(num_var, alpha, beta_, sigma_, p, T):
        (* Black Box call *)
        Y := [seq(B(Psi_alpha[i], p), i=1..T)]:
        M := Matrix(Y);
        row, col := Dimension(M):
        if row = 1 then
            lin_sys := false;
            t_helper := time():
            u := cppNewtonInterp(alpha,Y,x,p);
            t_cppNewton_total := t_cppNewton_total + (time() - t_helper):
            t_helper := time():
            result := [[MQRFR(m, u, 0, 1, p)]];
            t_mqrfr_total := t_mqrfr_total + (time() - t_helper):
            dq := result[1][3];
        else
            lin_sys := true;
            t_helper := time():
            u := get_u(M, col, alpha, p);
            t_cppNewton_total := t_cppNewton_total + (time() - t_helper):
        end if;
        if lin_sys then
            t_helper := time():
            for i from 1 to nops(u) do
                if MQRFR_done[i] then next; end if;
                temp := [op(temp), MQRFR(m, u[i], 0, 1, p)];
                result[i] := temp;
                temp := [];
            end do;
            t_mqrfr_total := t_mqrfr_total + (time() - t_helper):
            DQ := [seq(result[i][3], i=1..nops(result))];
            for i from 1 to numelems(DQ) do
                if DQ[i] > 1 then MQRFR_done[i] := true; end if;
            end do;
            dq := min(op(DQ));
        end if;
        if dq > 1 then
            print("NDSA: termination condition met");
            return result, lin_sys;
        else
            print("NDSA: MQRFR failed; doubling T");
            T := T*2;
            for i from 1 to num_eqn do
                if not MQRFR_done[i] then result[i] := []; end if;
            end do;
            DQ := [];
        end if;
        if T > 2^10 then
            print("NDSA: hit safety break at T=", T);
            return result, lin_sys;
        end if;
    end do;
end proc:

Deterministic_NDSA := proc(B, sigma_, beta_, num_var, p, num_points, max_points, num_eqn)
    local T, alpha, m, Psi_alpha, Y, u, dq, i, r,
          lin_sys, result, M, row, col;
    #print("In Deterministic_NDSA");
    lin_sys := false;
    T := max_points;
    result := [seq([], i=1..num_eqn)]:
    r := rand(p);
    alpha := [seq(i, i=2..T+1)]:
    Psi_alpha := get_point_on_affine_line(num_var, alpha, beta_, sigma_, p, T):
    Y := [seq(B(Psi_alpha[i], p), i=1..T)]:
    for i from 1 to num_eqn do
        m := Expand(product(x - alpha[j], j=1..num_points[i])) mod p:
        M := convert(Y[1..num_points[i]], Matrix);
        row, col := Dimension(M):
        if row = 1 then
            #u := Interp(alpha[1..num_points[i]], Y[1..num_points[i]], x) mod p:
            u := cppNewtonInterp(alpha[1..num_points[i]],Y[1..num_points[i]],x,p):
	else
            lin_sys := true;
            u := Deterministic_get_u(M, i, alpha[1..num_points[i]], p);
        end if;
        if not lin_sys then
            result := [MQRFR(m, u, 0, 1, p)];
            dq := result[1][3];
        else
            result[i] := [MQRFR(m, u, 0, 1, p)];
        end if;
    end do;
    return result;
end proc:

MRFI := proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, s, Primes, direction, sigma_, num_eval, den_eval, u, mon,
          numerator_done, denominator_done, Tcur, jDone,
          mqrfr_results, lin_sys, num_points_mqrfr,
          Numerators, Denominiators, deg_num, deg_den,
          lambda_num, lambda_den, terms_num, terms_den,
          R_num, R_den, Roots_num_eval, Roots_den_eval,
          num_mono, den_mono, coeff_num, coeff_den,
          final_num, final_den, temp, common_den_flag,
          bmea_done, temp_den, all_den_done, all_done, max_num_points_mqrfr,
          init_sigma,
          sampleCounts, mMax, tempG, ratReconVal,
          sigma_j, alphaVal, modulusTable, modulusPoly,
          Psi_alpha, BBvals, m_i, Y, interpVal, rr, tmpNum, tmpDen,
          t_helper;
    global mrfi_stats, bb_phase,
           t_cppNewton_total, t_cppRR_total, t_bmea_total, t_zippel_total;

    #print("MRFI -- num_vars=", num_vars, " num_eqn=", num_eqn, " p=", p);

    Primes := [seq(ithprime(i), i=1..num_vars)]:
    direction := [seq(i + 6, i=1..num_vars - 1)]:
    #print("MRFI direction = ", direction);

    (* Originally, the initial sigma array was [1,1,...,1] which makes a symmetric Toeplitz matrix
    degenerate to the all-ones matrix (rank 1, singular) so now we are using distinct
    primes so the first evaluation point is non-degenerate. *)

    init_sigma := Primes:

    sigma_ := []:
    num_eval := [seq([], i=1..num_eqn)]:
    den_eval := [seq([], i=1..num_eqn)]:
    numerator_done   := [seq(false, i=1..num_eqn)]:
    denominator_done := [seq(false, i=1..num_eqn)]:
    bmea_done := [seq(false, i=1..num_eqn)]:
    all_done := true:

    lambda_num := table(): terms_num := table(): R_num := table():
    lambda_den := table(): terms_den := table(): R_den := table():
    final_num  := table(): final_den := table():
    num_mono   := table(): coeff_num := table():
    den_mono   := table(): coeff_den := table():

    for i from 1 to num_eqn do
        lambda_num[i] := []: terms_num[i] := []: R_num[i] := []:
        lambda_den[i] := []: terms_den[i] := []: R_den[i] := []:
    end do:

    Tcur := 4:
    num_points_mqrfr := [seq(0, i=1..num_eqn)]:

    bb_phase := "NDSA":
    mqrfr_results, lin_sys := NDSA(B, init_sigma, direction,
                                    num_vars, p, Tcur, num_eqn):
    bb_phase := "MRFI":

    Numerators    := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    Denominiators := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:

    for k from 1 to num_eqn do
        num_eval[k] := []:
        den_eval[k] := []:
    end do:

    deg_num := [seq(degree(Numerators[i],    x), i=1..nops(Numerators))]:
    deg_den := [seq(degree(Denominiators[i], x), i=1..nops(Denominiators))]:
    for i from 1 to numelems(deg_den) do
        num_points_mqrfr[i] := deg_num[i] + deg_den[i] + 1:
    end do:
    max_num_points_mqrfr := max(op(deg_num)) + max(op(deg_den)) + 1:
    #print("MRFI num_points_mqrfr = ", num_points_mqrfr);
    #print("MRFI max_num_points_mqrfr = ", max_num_points_mqrfr);

    sampleCounts := num_points_mqrfr:
    mMax         := max_num_points_mqrfr:
    tempG        := rand(1..p-1):
    ratReconVal  := table():
    for i from 1 to num_eqn do ratReconVal[i] := table() end do:
    jDone := 0:

    common_den_flag := true:
    for k from 2 to num_eqn do
        if Denominiators[k] <> Denominiators[1] then
            common_den_flag := false; break;
        end if;
    end do;
    #print("MRFI common_den_flag = ", common_den_flag);

    while true do
        for j from jDone+1 to 2*Tcur do
            sigma_j := [seq(Primes[k]^j mod p, k=1..nops(Primes))]:
            sigma_  := [op(sigma_), sigma_j]:

            if lin_sys then
                alphaVal := [seq(tempG(), s=1..mMax)]:

                modulusTable := table():
                modulusPoly  := 1:
                for s from 1 to mMax do
                    modulusPoly       := Expand(modulusPoly*(x-alphaVal[s])) mod p:
                    modulusTable[s]   := modulusPoly:
                end do:
                Psi_alpha := get_point_on_affine_line(num_vars, alphaVal,
                                                     direction, sigma_j,
                                                     p, mMax):
                (* Black Box call *)
                BBvals := [seq(B(Psi_alpha[s], p), s=1..mMax)]:
                for i from 1 to num_eqn do
                    if numerator_done[i] and denominator_done[i] then next end if:
                    m_i := sampleCounts[i]:

                    Y         := [seq(BBvals[s][i], s=1..m_i)]:
                    t_helper := time():
                    interpVal := cppNewtonInterp([op(1..m_i, alphaVal)], Y, x, p):
                    t_cppNewton_total := t_cppNewton_total + (time() - t_helper):
                    t_helper := time():
                    rr        := cppRR(interpVal, modulusTable[m_i],
                                       x, deg_num[i], deg_den[i], p):
                    t_cppRR_total := t_cppRR_total + (time() - t_helper):

                    ratReconVal[i][j] := rr:

                    tmpNum := Eval(numer(rr), x = sigma_j[1]) mod p:
                    tmpDen := Eval(denom(rr), x = sigma_j[1]) mod p:

                    if not numerator_done[i] then
                        num_eval[i] := [op(num_eval[i]), tmpNum]:
                    end if:
                    if not denominator_done[i] then
                        den_eval[i] := [op(den_eval[i]), tmpDen]:
                    end if:
                end do:
            else
                bb_phase := "NDSA":
                mqrfr_results, lin_sys := NDSA(B, sigma_j, direction,
                                                num_vars, p,
                                                max_num_points_mqrfr, 1):
                bb_phase := "MRFI":
                Numerators    := [seq(mqrfr_results[k][1],
                                      k=1..nops(mqrfr_results))]:
                Denominiators := [seq(mqrfr_results[k][2],
                                      k=1..nops(mqrfr_results))]:
                for k from 1 to num_eqn do
                    if not numerator_done[k] then
                        num_eval[k] := [op(num_eval[k]),
                                        eval(Numerators[k], x=sigma_j[1]) mod p]:
                    end if;
                    if not denominator_done[k] then
                        den_eval[k] := [op(den_eval[k]),
                                        eval(Denominiators[k], x=sigma_j[1]) mod p]:
                    end if;
                end do;
            end if;
        end do;

        jDone := 2*Tcur:

        all_done := true:
        for k from 1 to num_eqn do
            if numerator_done[k] then next; end if;
            t_helper := time():
            lambda_num[k] := BMEA(num_eval[k], p, Z):
            t_bmea_total := t_bmea_total + (time() - t_helper):
            terms_num[k]  := degree(lambda_num[k], Z):
            R_num[k]      := Roots(lambda_num[k]) mod p:
            if R_num[k] = [] then numerator_done[k] := false; next; end if;
            if nops(R_num[k]) > 0 and R_num[k][1][1] = 0 then
                R_num[k] := remove(x -> x = [0,1], R_num[k]):
                terms_num[k] := terms_num[k] - 1:
            end if;
            if nops(R_num[k]) = terms_num[k]
               and terms_num[k] < iquo(nops(num_eval[k]), 2) then
                numerator_done[k] := true:
            end if;
        end do;

        all_den_done := true:
        if common_den_flag then
            if not denominator_done[1] then
                t_helper := time():
                lambda_den[1] := BMEA(den_eval[1], p, Z):
                t_bmea_total := t_bmea_total + (time() - t_helper):
                terms_den[1]  := degree(lambda_den[1], Z):
                R_den[1]      := Roots(lambda_den[1]) mod p:
                if nops(R_den[1]) > 0 and R_den[1][1][1] = 0 then
                    R_den[1] := remove(x -> x = [0,1], R_den[1]):
                    terms_den[1] := terms_den[1] - 1:
                end if;
                if nops(R_den[1]) = terms_den[1]
                   and terms_den[1] < iquo(nops(den_eval[1]), 2) then
                    for k from 1 to num_eqn do
                        denominator_done[k] := true:
                    end do;
                end if;
            end if;
        else
            for k from 1 to num_eqn do
                if denominator_done[k] then next; end if;
                t_helper := time():
                lambda_den[k] := BMEA(den_eval[k], p, Z):
                t_bmea_total := t_bmea_total + (time() - t_helper):
                terms_den[k]  := degree(lambda_den[k], Z):
                R_den[k]      := Roots(lambda_den[k]) mod p:
                if nops(R_den[k]) > 0 and R_den[k][1][1] = 0 then
                    R_den[k] := remove(x -> x = [0,1], R_den[k]):
                    terms_den[k] := terms_den[k] - 1:
                end if;
                if nops(R_den[k]) = terms_den[k]
                   and terms_den[k] < iquo(nops(den_eval[k]), 2) then
                    denominator_done[k] := true:
                end if;
            end do;
        end if;

        for i from 1 to num_eqn do
            bmea_done[i] := numerator_done[i] and denominator_done[i]:
        end do;
        all_done := true;
        for i from 1 to num_eqn do all_done := all_done and bmea_done[i]; end do;
        #print("MRFI bmea_done = ", bmea_done, " all_done = ", all_done);
        if all_done then break; end if;

        Tcur := 2*Tcur:
        if Tcur > 2^14 then
            print("MRFI: hit safety break at Tcur=", Tcur);
            break;
        end if;
    end do;

    Roots_num_eval := [seq([seq(r[1], r in R_num[k])], k=1..num_eqn)]:
    if common_den_flag then
        Roots_den_eval := [[seq(r[1], r in R_den[1])]]:
        for k from 2 to num_eqn do
            Roots_den_eval := [op(Roots_den_eval), Roots_den_eval[1]]:
        end do;
    else
        Roots_den_eval := [seq([seq(r[1], r in R_den[k])], k=1..num_eqn)]:
    end if;

    for k from 1 to num_eqn do
        temp := generate_monomials(Roots_num_eval[k], num_vars, Primes, vars):
        if temp = FAIL then return FAIL; end if;
        num_mono[k] := temp;
        t_helper := time():
        coeff_num[k] := Zippel_Transpose_Vandermonde_solver(num_eval[k],
                              terms_num[k], Roots_num_eval[k], lambda_num[k], p):
        t_zippel_total := t_zippel_total + (time() - t_helper):
    end do;

    if common_den_flag then
        temp := generate_monomials(Roots_den_eval[1], num_vars, Primes, vars):
        if temp = FAIL then return FAIL; end if;
        den_mono[1]  := temp;
        t_helper := time():
        coeff_den[1] := Zippel_Transpose_Vandermonde_solver(den_eval[1],
                              terms_den[1], Roots_den_eval[1], lambda_den[1], p):
        t_zippel_total := t_zippel_total + (time() - t_helper):
        for k from 2 to num_eqn do
            den_mono[k]  := den_mono[1]:
            coeff_den[k] := coeff_den[1]:
        end do;
    else
        for k from 1 to num_eqn do
            temp := generate_monomials(Roots_den_eval[k], num_vars, Primes, vars):
            if temp = FAIL then return FAIL; end if;
            den_mono[k] := temp;
            t_helper := time():
            coeff_den[k] := Zippel_Transpose_Vandermonde_solver(den_eval[k],
                              terms_den[k], Roots_den_eval[k], lambda_den[k], p):
            t_zippel_total := t_zippel_total + (time() - t_helper):
        end do;
    end if;

    for k from 1 to num_eqn do
        lcoeff(add(mon, mon in den_mono[k]), vars, 'mon');
        if not member(mon, den_mono[k], 'i') then
            error "bug in leading monomial";
        end if;
        u := 1/coeff_den[k][i] mod p;
        coeff_num[k] := u * coeff_num[k] mod p:
        coeff_den[k] := u * coeff_den[k] mod p:
        final_num[k] := construct_final_polynomial(coeff_num[k], num_mono[k]):
        final_den[k] := construct_final_polynomial(coeff_den[k], den_mono[k]):
    end do;

    mrfi_stats := table():
    mrfi_stats["terms_num"]       := [seq(terms_num[k], k=1..num_eqn)]:
    mrfi_stats["terms_den"]       := [seq(terms_den[k], k=1..num_eqn)]:
    mrfi_stats["deg_num"]         := deg_num:
    mrfi_stats["deg_den"]         := deg_den:
    mrfi_stats["common_den_flag"] := common_den_flag:
    mrfi_stats["mMax"]            := mMax:
    mrfi_stats["sampleCounts"]    := sampleCounts:

    return final_num, final_den;
end proc:

RandRational := proc(N::posint)
    return proc() local a, b;
        a := rand(-N..N)();
        b := rand(1..N)();
        if a = 0 then 0 else a/b fi;
    end proc;
end proc:

get_data := proc(test_case)
    local Sys, Vars, i, ii, jj, params, ff, gg, n;
    if nargs = 1 then
        if test_case = "bspline" then
            Sys := {x7 + x12 - 1, x8 + x13 - 1, x21 + x6 + x11 - 1,
                    x1*y1 + x1 - x2, x11*y3 + x11 - x12, x16*y5 - x17*y5 - x17,
                    -x20*y3 + x21*y3 + x21, x3*y2 + x3 - x4,
                    -x8*y4 + x9*y3 + x9, 2*x1*y1^2 - 2*x1 - 2*x10 + 4*x2,
                    -x10*y2 + x18*y2 + x18 - x19, 2*x11*y3^2 - 2*x11 + 4*x12 - 2*x13,
                    -x13*y4 + x14*y4 + x14 - x15, 2*x15*y5^2 - 4*x16*y5^2 + 2*x17*y5^2 - 2*x17,
                    2*x19*y3^2 - 4*x20*y3^2 + 2*x21*y3^2 - 2*x21,
                    2*x3*y2^2 - 2*x3 + 4*x4 - 2*x5, -x5*y3 + x6*y3 + x6 - x7,
                    2*x7*y4^2 - 4*x8*y4^2 + 2*x9*y4^2 - 2*x9,
                    -4*x10*y2^2 + 2*x18*y2^2 + 2*x2*y2^2 - 2*x18 + 4*x19 - 2*x20,
                    2*x12*y4^2 - 4*x13*y4^2 + 2*x14*y4^2 - 2*x14 + 4*x15 - 2*x16,
                    2*x4*y3^2 - 4*x5*y3^2 + 2*x6*y3^2 - 2*x6 + 4*x7 - 2*x8}:
        elif test_case = "small_sys_low_deg" then
            Sys := {x1+y1*x2+y2*x3-1, y2*x1+x2+y1*x3-2, (y1-y2)*x1-x2+y2*x3-7}:
        elif test_case = "small_Sys" then
            Sys := {x1+y1*x2+y1-3, y2*x1+x2+y1-1}:
        elif test_case = "mike" then
            Sys := {y1*x1+y1*x2-1, y1*y2*x1-x2-1}:
        elif test_case = "example" then
            Sys := {(y1*y2-1)*x1 + (y1^2-2*y1+3),
                    (y1*y2-1)*x2 + (y1*y2-y1-3*y2+1)}:
        elif test_case = "bsbug" then
            Sys := {(2*y3^2*y4 - y3*y4^2 + 3*y3*y4 - y4^2 + y3 + y4 + 1)*x1
                    = y3*y4^2}:
        elif test_case = 1 then
            ff := y; gg := x - 4; Vars := [x, y]:
            return Vars, ff, gg, numelems(Vars), 1, Vars;
        end if;
    elif nargs > 1 then
        Vars := [seq(x||i, i=1..args[2])]:
        if test_case = "rand" then
            ff := randpoly(Vars, terms=args[3]):
            gg := randpoly(Vars, terms=args[4]):
            return Vars, ff, gg, numelems(Vars), 1, Vars;
        elif test_case = "rat_rand" then
            ff := randpoly(Vars, coeffs=RandRational(args[5]), terms=args[3]):
            gg := randpoly(Vars, coeffs=RandRational(args[6]), terms=args[4]):
            return Vars, ff, gg, numelems(Vars), 1, Vars;
        elif test_case = "TP" then
            n := args[2]:
            Vars   := [seq(x||i, i=1..n)]:
            params := [seq(y||i, i=1..n)]:
            Sys := { seq( add( y||(abs(ii-jj)+1) * x||jj, jj=1..n ) - 1,
                          ii = 1..n ) }:
            return Sys, Vars, params, numelems(params), numelems(Vars);
        end if;
    end if;

    # Fallback for the original system test cases (mike, bspline, etc.)
    # F1: Vars now returned as an ordered LIST and params as an ordered LIST.
    Vars := [seq(x||i, i=1..nops(Sys))]:
    params := convert(indets(Sys) minus convert(Vars, set), list):
    return Sys, Vars, params, nops(params), nops(Vars);
end proc:

test_prime := 2^31 - 1:  
n_min := 4:
n_max := 12:
do_verify := true:
do_ffge := true:
summary := []:

for n_test from n_min to n_max do
    #print("=================================================================="):
    #printf("  TEST: symmetric Toeplitz, n = %d\n", n_test):
    #print("=================================================================="):

    Sys, Vars, params, num_vars, num_eqn := get_data("TP",n_test):
    #print("Vars   = ", Vars):
    #print("params = ", params):
    #print("|Sys|  = ", nops(Sys)):

    counter         := 0:
    num_lines       := 0:
    bb_phase        := "MRFI":       # default phase; flipped to "NDSA" around NDSA calls
    bb_bench_calls  := 1000:         # iterations per BB call for inline timing
    bb_calls_ndsa   := 0:
    bb_calls_mrfi   := 0:
    t_ndsa_total    := 0.0:
    t_mrfi_total    := 0.0:
    t_cppNewton_total := 0.0:
    t_mqrfr_total     := 0.0:
    t_cppRR_total     := 0.0:
    t_bmea_total      := 0.0:
    t_zippel_total    := 0.0:
    t_iratrecon_total := 0.0:
    forget(get_point_on_affine_line): 

    B := Constuct_Sys_Blackbox(Sys, Vars, params):

    Num, Den := 'Num', 'Den':
    status := "OK":
    mrfi_out := []:
    t_mrfi_wall_start := time():
    try
        mrfi_out := [MRFI(B, num_vars, num_eqn, params, test_prime)]:
    catch:
        status := cat("ERROR: ", StringTools:-FormatMessage(lastexception[2..-1])):
        print("MRFI threw: ", status):
    end try;
    t_mrfi_wall := time() - t_mrfi_wall_start:
    mrfi_calls := counter:
    # bb_calls_ndsa, bb_calls_mrfi, t_ndsa_total, t_mrfi_total already populated
    if status = "OK" then
        if nops(mrfi_out) = 2 and mrfi_out <> [FAIL, FAIL] then
            Num := mrfi_out[1]:
            Den := mrfi_out[2]:
        else
            status := "MRFI returned FAIL":
        end if;
    end if;

    if status = "OK" then
        stats_terms_num := mrfi_stats["terms_num"]:
        stats_terms_den := mrfi_stats["terms_den"]:
        stats_deg_num   := mrfi_stats["deg_num"]:
        stats_deg_den   := mrfi_stats["deg_den"]:
        stats_common    := mrfi_stats["common_den_flag"]:
        stats_mMax      := mrfi_stats["mMax"]:
    else
        stats_terms_num := []: stats_terms_den := []:
        stats_deg_num   := []: stats_deg_den   := []:
        stats_common    := false: stats_mMax := 0:
    end if;

    if do_ffge then
        Y_ffge := [seq(y||i, i=1..n_test)]:
        A_ffge := LinearAlgebra:-ToeplitzMatrix(Y_ffge, symmetric):
        b_ffge := Vector(n_test, fill=1):
        ffge_status := "OK":
        try
            det_ffge, f_ffge, g_ffge := FFGE(A_ffge, b_ffge, Y_ffge):
        catch:
            ffge_status := cat("ERROR: ",
                               StringTools:-FormatMessage(lastexception[2..-1])):
            print("FFGE threw: ", ffge_status):
        end try:

        if ffge_status = "OK" then
            ffge_terms_num   := ffge_stats["f_terms"]:
            ffge_terms_den   := ffge_stats["g_terms"]:
            ffge_y_terms     := ffge_stats["y_terms"]:
            ffge_det_terms   := ffge_stats["det_terms"]:
            ffge_elim_swell  := ffge_stats["elim_num_max"]:
            ffge_backsub_swell := ffge_stats["backsub_N_max"]:

            # FFGE timing benchmark — same methodology as BB above.
            # Run FFGE ffge_bench_calls times at the same input, average.
            ffge_bench_calls := 1:
            t_ffge_start := time():
            to ffge_bench_calls do
                FFGE(A_ffge, b_ffge, Y_ffge):
            end do:
            t_ffge := evalf((time() - t_ffge_start) / ffge_bench_calls):
        else
            ffge_terms_num     := []:  ffge_terms_den     := []:
            ffge_y_terms       := []:  ffge_det_terms     := -1:
            ffge_elim_swell    := -1:  ffge_backsub_swell := -1:
            t_ffge             := 0.0:
        end if:
    else
        ffge_status        := "SKIPPED":
        t_ffge             := 0.0:
        ffge_terms_num     := []:  ffge_terms_den     := []:
        ffge_y_terms       := []:  ffge_det_terms     := -1:
        ffge_elim_swell    := -1:  ffge_backsub_swell := -1:
    end if;

    if status = "OK" then
        Ratrecon_num   := table():
        Ratrecon_den   := table():
        Final_rat_poly := table():
        t_helper := time():
        for i from 1 to num_eqn do
            Ratrecon_num[i] := iratrecon(Num[i], test_prime):
            Ratrecon_den[i] := iratrecon(Den[i], test_prime):
            Final_rat_poly[i] := Ratrecon_num[i] / Ratrecon_den[i]:
        end do;
        t_iratrecon_total := time() - t_helper:

        if do_verify then
            try
                A_ref, b_ref := GenerateMatrix(Sys, Vars):
                x_ref := LinearSolve(A_ref, b_ref):
                all_match := true:
                for i from 1 to num_eqn do
                    diff_i := normal(Final_rat_poly[i] - x_ref[i]):
                    if diff_i <> 0 then
                        all_match := false:
                        #printf("  x%d MISMATCH: residue = %a\n",
                               #i, diff_i):
                        #printf("    recovered = %a\n", Final_rat_poly[i]):
                        #printf("    reference = %a\n", x_ref[i]):
                    end if;
                end do;
                if all_match then
                    printf("  PASS  n=%2d   "
                           "BB(NDSA)=%d  BB(MRFI)=%d  "
                           "t_NDSA=%.3f s  t_MRFI=%.3f s  overall=%.3f s\n",
                           n_test, bb_calls_ndsa, bb_calls_mrfi,
                           t_ndsa_total, t_mrfi_total,
                           t_ndsa_total + t_mrfi_total):
                    summary := [op(summary),
                                [n_test, "PASS", mrfi_calls,
                                 bb_calls_ndsa,   bb_calls_mrfi,
                                 stats_terms_num, stats_terms_den,
                                 stats_deg_num,   stats_deg_den,
                                 stats_common,    stats_mMax,
                                 ffge_status,     t_ffge,
                                 ffge_terms_num,  ffge_terms_den,
                                 ffge_y_terms,    ffge_det_terms,
                                 ffge_elim_swell, ffge_backsub_swell,
                                 t_ndsa_total,    t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total, t_mqrfr_total, t_cppRR_total,
                                 t_bmea_total,      t_zippel_total, t_iratrecon_total]]:
                else
                    printf("  FAIL  n=%2d (mismatch)\n", n_test):
                    summary := [op(summary),
                                [n_test, "FAIL-mismatch", mrfi_calls,
                                 bb_calls_ndsa,   bb_calls_mrfi,
                                 stats_terms_num, stats_terms_den,
                                 stats_deg_num,   stats_deg_den,
                                 stats_common,    stats_mMax,
                                 ffge_status,     t_ffge,
                                 ffge_terms_num,  ffge_terms_den,
                                 ffge_y_terms,    ffge_det_terms,
                                 ffge_elim_swell, ffge_backsub_swell,
                                 t_ndsa_total,    t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total, t_mqrfr_total, t_cppRR_total,
                                 t_bmea_total,      t_zippel_total, t_iratrecon_total]]:
                end if;
            catch:
                printf("  PARTIAL  n=%2d   (ref solve threw)\n", n_test):
                for i from 1 to num_eqn do
                    printf("    x%d = %a\n", i, Final_rat_poly[i]):
                end do;
                summary := [op(summary),
                            [n_test, "NO-VERIFY", mrfi_calls,
                             bb_calls_ndsa,   bb_calls_mrfi,
                             stats_terms_num, stats_terms_den,
                             stats_deg_num,   stats_deg_den,
                             stats_common,    stats_mMax,
                                 ffge_status,     t_ffge,
                                 ffge_terms_num,  ffge_terms_den,
                                 ffge_y_terms,    ffge_det_terms,
                                 ffge_elim_swell, ffge_backsub_swell,
                                 t_ndsa_total,    t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total, t_mqrfr_total, t_cppRR_total,
                                 t_bmea_total,      t_zippel_total, t_iratrecon_total]]:
            end try;
        else
            printf("  RECOVERED (unverified) n=%2d  BB-calls=%d\n",
                   n_test, mrfi_calls):
            for i from 1 to num_eqn do
                printf("    x%d = %a\n", i, Final_rat_poly[i]):
            end do;
            summary := [op(summary),
                        [n_test, "UNVERIFIED", mrfi_calls,
                         bb_calls_ndsa,   bb_calls_mrfi,
                         stats_terms_num, stats_terms_den,
                         stats_deg_num,   stats_deg_den,
                         stats_common,    stats_mMax,
                                 ffge_status,     t_ffge,
                                 ffge_terms_num,  ffge_terms_den,
                                 ffge_y_terms,    ffge_det_terms,
                                 ffge_elim_swell, ffge_backsub_swell,
                                 t_ndsa_total,    t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total, t_mqrfr_total, t_cppRR_total,
                                 t_bmea_total,      t_zippel_total, t_iratrecon_total]]:
        end if;
    else
        printf("  FAIL  n=%2d   %s\n", n_test, status):
        summary := [op(summary),
                    [n_test, status, mrfi_calls,
                     bb_calls_ndsa,   bb_calls_mrfi,
                     stats_terms_num, stats_terms_den,
                     stats_deg_num,   stats_deg_den,
                     stats_common,    stats_mMax,
                                 ffge_status,     t_ffge,
                                 ffge_terms_num,  ffge_terms_den,
                                 ffge_y_terms,    ffge_det_terms,
                                 ffge_elim_swell, ffge_backsub_swell,
                                 t_ndsa_total,    t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total, t_mqrfr_total, t_cppRR_total,
                                 t_bmea_total,      t_zippel_total, t_iratrecon_total]]:
    end if;
end do;

print("=================================================================="):
print("  SUMMARY"):
print("=================================================================="):
printf("%4s  %-12s  %8s  %8s  %10s  %10s  %14s  %14s  %14s\n",
       "n", "status",
       "BB-NDSA", "BB-MRFI",
       "t-NDSA(s)", "t-MRFI(s)",
       "wall-raw(s)", "wall-corr(s)", "FFGE(s)"):
for entry in summary do
    printf("%4d  %-12s  %8d  %8d  %10.3f  %10.3f  %14.3f  %14.3f  %14.6f\n",
           entry[1], entry[2],
           entry[4], entry[5],
           entry[20], entry[21],
           entry[22], entry[22] - 999.0 * (entry[20] + entry[21]),
           entry[13]):
end do;

print(""):
print("  Term count comparison: MRFI vs Lipson FFGE"):
printf("%4s  %-22s  %-22s  %-22s  %-22s  %10s\n",
       "n",
       "MRFI num terms", "FFGE f terms",
       "MRFI den terms", "FFGE g terms",
       "FFGE max swell"):
for entry in summary do
    printf("%4d  %-22a  %-22a  %-22a  %-22a  %10a\n",
           entry[1], entry[6], entry[14], entry[7], entry[15],
           max(entry[18], entry[19])):
end do;

report_path := "TimingsTN.txt":
fd := fopen(report_path, WRITE):
fprintf(fd, "============================================================\n"):
fprintf(fd, "  Symmetric Toeplitz MRFI benchmark\n"):
fprintf(fd, "  Prime p = %d\n", test_prime):
fprintf(fd, "  Range   n = %d .. %d\n", n_min, n_max):
fprintf(fd, "============================================================\n\n"):

for entry in summary do
    fprintf(fd, "  n = %d\n",                  entry[1]):
    fprintf(fd, "  MRFI Status                           : %s\n",     entry[2]):
    fprintf(fd, "  Total BB calls                        : %d  (NDSA = %d, MRFI = %d)\n",
            entry[3], entry[4], entry[5]):
    fprintf(fd, "  BB calls in NDSA                      : %d\n",     entry[4]):
    fprintf(fd, "  BB calls in MRFI                      : %d\n",     entry[5]):
    fprintf(fd, "  Total NDSA time for BB calls (s)      : %.9f\n",
            entry[20]):
    fprintf(fd, "  Total MRFI time for BB calls (s)      : %.9f\n",   entry[21]):
    fprintf(fd, "  Total BB calls time (NDSA+MRFI)       : %.9f\n",
            entry[20] + entry[21]):
    fprintf(fd, "  Total time (s)                        : %.9f \n",
            entry[22] - 999.0 * (entry[20] + entry[21])):
    fprintf(fd, "  Time for MRFI Individual component - \n"):
    fprintf(fd, "  Time -> cppNewtonInterp (s)  : %.9f\n",   entry[23]):
    fprintf(fd, "  Time -> MQRFR (s)            : %.9f\n",   entry[24]):
    fprintf(fd, "  Time -> cppRR (s)            : %.9f\n",   entry[25]):
    fprintf(fd, "  Time -> BMEA (s)             : %.9f\n",   entry[26]):
    fprintf(fd, "  Time -> Vandermonde (s)      : %.9f\n",   entry[27]):
    fprintf(fd, "  Deg_num (per equation)   : %a\n",     entry[8]):
    fprintf(fd, "  Deg_den (per equation)   : %a\n",     entry[9]):
    fprintf(fd, "  Terms_num (per equation) : %a\n",     entry[6]):
    fprintf(fd, "  Terms_den (per equation) : %a\n",     entry[7]):
    fprintf(fd, " \n"):
    fprintf(fd, "  FFGE status              : %s\n",     entry[12]):
    fprintf(fd, "  FFGE total time (s)      : %.9f\n",
            entry[13]):
    fprintf(fd, "  FFGE f[i] terms          : %a\n",     entry[14]):
    fprintf(fd, "  FFGE g[i] terms          : %a\n",     entry[15]):
    fprintf(fd, "  FFGE y[i] terms (Pre GCD): %a\n",     entry[16]):
    fprintf(fd, "  FFGE det(A) terms        : %a\n",     entry[17]):
    fprintf(fd, "  FFGE Max Elim. step swell: %a  (numterms(num) before exact div)\n",
            entry[18]):
    fprintf(fd, "  FFGE max Back sub swell  : %a  (numterms(N[i]) before exact div)\n",
            entry[19]):
    fprintf(fd, "  Term count comparison (MRFI vs FFGE)\n"):
    fprintf(fd, "  Numerator terms MRFI     : %a\n",     entry[6]):
    fprintf(fd, "  Numerator terms FFGE     : %a\n",     entry[14]):
    fprintf(fd, "  Denominator terms MRFI   : %a\n",     entry[7]):
    fprintf(fd, "  Denominator terms FFGE   : %a\n\n",   entry[15]):
end do:
fclose(fd):
