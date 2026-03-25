construct_final_polynomial:=proc(coeff_,Monomials)
    print("In construct_final_polynomial"):
    local i,f:
    f:=0:
    for i from 1 to nops(coeff_) do
        f:=f+coeff_[i]*Monomials[i]:
         print("YES:"):
    end do:
    return f:
end proc:
