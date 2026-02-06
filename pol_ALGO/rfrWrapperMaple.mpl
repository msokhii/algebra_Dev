lib := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/wrapOBJ.so":

mRATRECON := define_external(
    'ratRECON_C',
    mLen::integer[4],
    degM::integer[4],
    M::ARRAY(0..mLen-1, datatype=integer[8]),
    uLen::integer[4],
    degU::integer[4],
    U::ARRAY(0..uLen-1, datatype=integer[8]),
    N::integer[4],
    DBound::integer[4],
    p::integer[8],
    nOLEN::integer[4],
    nOUT::ARRAY(0..nOLEN-1, datatype=integer[8]),
    degNOUT::REF(integer[4]),
    dOLEN::integer[4],
    dOUT::ARRAY(0..dOLEN-1, datatype=integer[8]),
    degDOUT::REF(integer[4]),
    RETURN::integer[4],
    LIB=lib
):


# EXAMPLE: 

(* degU := 4: degM := 5:
mLen := degM+1: uLen := degU+1:

M := Array(0..mLen-1, datatype=integer[8], fill=0):
U := Array(0..uLen-1, datatype=integer[8], fill=0):

U[0]:=3: U[1]:=5: U[2]:=7: U[3]:=11: U[4]:=13:
M[5]:=1:

p := 97:
N := 2: DBound := 2:

# OUTPUT BUFFERS: 

nOLEN := N+1:
dOLEN := DBound+1:
nOUT := Array(0..nOLEN-1, datatype=integer[8], fill=0):
dOUT := Array(0..dOLEN-1, datatype=integer[8], fill=0):
degNOUT := 0:
degDOUT := 0:

rc := mRATRECON(mLen,degM,M,uLen,degU,U,N,DBound,p,nOLEN,nOUT,'degNOUT',dOLEN,dOUT,'degDOUT'):
degNOUT;
degDOUT;

print("rc", rc):
print("nOUT full", [seq(nOUT[i], i=0..nOLEN-1)]):
print("dOUT full", [seq(dOUT[i], i=0..dOLEN-1)]):
M;
U; *)

# TESTING PROCEDURE 1: 
# Something wrong with deg>60.

d := 20000;
pp := 2305843009213693951;
printf("DEGREE CHOSEN:%d\n",d);
NN := randpoly(x,degree=d,dense) mod pp;
DD := randpoly(x,degree=d,dense) mod pp:
while Gcd(NN,DD) mod pp <> 1 do 
    DD := randpoly(x,degree=d,dense) mod pp:
od;
DD := DD/lcoeff(DD) mod pp;
MM := mul((x-i),i=1..2*d+1);
MM := Expand(MM) mod pp:
g := Gcdex(DD,MM,x,'SS') mod pp:
printf("gcd(DD,MM) mod 101 = %a\n", g);
UU := Rem(SS*NN,MM,x) mod pp:
degU := degree(UU);
degM := degree(MM);
printf("TRUE DEGREE OF U:%d\nDEGREE OF M:%d\n",degU,degM);
M := Array(0..degM, [seq(coeff(MM,x,i),i=0..degM)], datatype=integer[8]):
U := Array(0..degU, [seq(coeff(UU,x,i),i=0..degU)], datatype=integer[8]): 
nOUT := Array(0..degM, datatype=integer[8] ):
dOUT := Array(0..degM, datatype=integer[8] ):
rc := mRATRECON(degM+1,degM,M,degU+1,degU,U,d,d,pp,d+1,nOUT,'degNOUT',d+1,dOUT,'degDOUT'):
degDOUT;
degNOUT;
printf("COMPUTED DEGREE OF N:%d\nDEGREE OF D:%d\n",degNOUT,degDOUT);
NNarr := Array(0..degM, [seq(coeff(NN, x, i), i=0..degNOUT)], datatype=integer[8]):
DDarr := Array(0..degM, [seq(coeff(DD, x, i), i=0..degDOUT)], datatype=integer[8]):
nOUT-NNarr;
dOUT-DDarr;

