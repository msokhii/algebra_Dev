get_point_on_affine_line_2 := proc(num_var, alpha, beta_, sigma_pt, p, T)
    local psi, np, nv;

    psi := Array(1..T);

    for np from 1 to T do
        psi[np] := Array(1..num_var);

        # First coordinate is the parameter value
        psi[np][1] := alpha[np];

        # Remaining coordinates
        for nv from 2 to num_var do
            psi[np][nv] := (beta_[nv-1]*alpha[np]
                          - beta_[nv-1]*sigma_pt[1]
                          + sigma_pt[nv]) mod p;
        end do;
    end do;

    return [seq(convert(psi[np], list), np=1..T)];
end proc: