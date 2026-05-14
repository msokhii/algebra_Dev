read "./mapleWrapperv2.mpl":

BMEA := proc(v::list, p::posint, Z::name)
    # Wraps the C++ Berlekamp-Massey via cppBMM, lifting the returned
    # coefficient list [Lambda_0, Lambda_1, ..., Lambda_deg] into a polynomial
    # in Z. Empty / all-zero input returns 1 so Roots(...) and degree(...)
    # downstream stay well-behaved.
    local L, d, i;
    if numelems(v) = 0 then return 1; end if;
    L := cppBMM(v, p);
    print(L);
    if L = [] then return 1; end if;
    d := nops(L) - 1;
    return add(L[i+1]*Z^i, i=0..d);
end:
# --- standalone BMEA sanity check ---
testBM := proc()
local p, F, v, BMEA_old, lam_old, lam_new;
    p := 2^31 - 1;
    # Fibonacci-like LFSR: a_{i+2} = a_{i+1} + a_i; Lambda(Z) = Z^2 - Z - 1
    v := [1, 1, 2, 3, 5, 8, 13, 21];

    # the ORIGINAL Maple BMEA, inlined here for comparison
    BMEA_old := proc(v, p, Z)
       local n, m, R0, R1, V0, V1, Q, ii;
       n  := iquo(nops(v), 2);
       m  := 2*n - 1;
       R0 := Z^(2*n);
       R1 := add(v[m+1-ii]*Z^ii, ii=0..m) mod p;
       V0 := 0; V1 := 1;
       while n <= degree(R1, Z) do
          R0, R1 := R1, Rem(R0, R1, Z, 'Q') mod p;
          V0, V1 := V1, Expand(V0 - Q*V1) mod p;
       end do;
       ii := 1/lcoeff(V1, Z) mod p;
       return ii*V1 mod p;
    end;

    lam_old := BMEA_old(v, p, Z);
    lam_new := BMEA(v, p, Z);                  # the new one from TMv1.mpl
    printf("old Lambda = %a\n",   lam_old);
    printf("new Lambda = %a\n",   lam_new);
    printf("difference  = %a\n",  Expand(lam_new - lam_old) mod p);
    printf("cppBMM raw  = %a\n",  cppBMM(v, p));
end:
testBM();
