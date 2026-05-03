read "./0_construct_BB.mpl":
read "./1_MRFI.mpl":

(* This is ver. 2 of MRFI2 - In use currently. *)
read "./MRFI2.mpl":
read "./2a_NDSA.mpl":
read "./2b_Deterministic_NDSA.mpl":
read "./3_get_point_on_affine_line.mpl":
read "./4_MQRFR.mpl":
read "./5_BMEA.mpl":
read "./6_generate_monomials.mpl":
read "./7_zippel_vandermonde_solver.mpl":
read "./8_construct_final_polynomial.mpl":
read "./9_rvr.mpl":
read "./data_gen.mpl":
read "./helpers.mpl":

(* CPP RATRECON WRAPPER OLD *)
# read "./cppRatRecon.mpl":
# read "./cppNewtonInterp.mpl":

(* CPP RATRECON + NEWTONINTERP WRAPPER OLD *)
read "./mapleWrapper.mpl":

(* SAMPLE TEST CASES *)
# test_case:="example":
# test_case:="small_sys_low_deg":
# test_case:="bspline":
#  test_case:="small_Sys":
# test_case:="mike":
# test_case:="bsbug":

matSize := 2:
maxTS := 14:
BBCalls := table():
for k from 1 to maxTS do
termsN := 0:
termsD := 0:
test_case := "TS":
num_lines:=0:
Sys, Vars, params, num_vars, num_eqn:= get_data(test_case,matSize):
counter := 0:
B := Constuct_Sys_Blackbox(Sys, Vars, params):
p:= 2^31-1:
try
Num,Den := rrMRFI(B, num_vars, num_eqn, params, p):
    catch:
    lprint("ERROR:", lasterror()):
end try:
print("NUM",Num);
print("DENUM",Den);
Ratrecon_num:=table():
Ratrecon_den:=table():
Final_rat_poly:=table():
for i from 1 to num_eqn do 
    Ratrecon_num[i]:=iratrecon(Num[i],p):
    #print("numerator = ",Ratrecon_num[i]):
    Ratrecon_den[i]:=iratrecon(Den[i],p):
    #print("denominator = ",Ratrecon_den[i]):
end do:

#print("IRAT NUM",convert(Ratrecon_num,list));
#print("IRAT DENUM",convert(Ratrecon_den,list));

print("======================================================"):
print("Displaying the results"):

if(num_eqn >1)then 
    og_soln:=get_eqn(Sys,Vars):
    og_soln:=convert(og_soln,list):
    for i from 1 to num_eqn do 
        print("x",i,"="):
        Final_rat_poly[i]:=Ratrecon_num[i]/Ratrecon_den[i]:
        termsN := termsN+nops(Ratrecon_num[i]):
        termsD := termsD+nops(Ratrecon_den[i]):
        lprint("Recovered Polynomial = ",Final_rat_poly[i]):
        lprint("Original Polynomial  = ",op(2,og_soln[i])):
        printf("f%d/g%d-ff%d/gg%d = %a\n",i,i,i,i,simplify(Final_rat_poly[i]-op(2,og_soln[i])));
    end do:
    elif num_eqn =1 then 
        Final_rat_poly[1]:=Ratrecon_num[1]/Ratrecon_den[1]:
        termsN := termsN+nops(Ratrecon_num[i]):
        termsD := termsD+nops(Ratrecon_den[i]):
        lprint("Recovered Polynomial = ",Final_rat_poly[1]):
        lprint("Original polynomial  = ",F/G):
        printf("f1/g1 - F/G = %a\n",simplify(Final_rat_poly[1]-F/G));
end if:

print("======================================================"):
print("Total number of lines generated in get_point_on_affine_line:", num_lines):
lprint("Total Black Box Calls:", counter):
BBCalls[k] := [matSize,counter,floor(evalf(termsN/num_eqn)),floor(evalf(termsD/num_eqn)),num_eqn,num_vars]:
matSize++:
od:

BBCalls := convert(BBCalls,list):
print(BBCalls);
