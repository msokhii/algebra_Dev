lib := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/wrapOBJ2.so":

mVSOLVE := define_external(
    'vSOLVER_C',
    RR::ARRAY(1..nn,datatype=integer[8]),
    LL::ARRAY(1..nn,datatype=integer[8]),
    nn::integer[4],
    aa::ARRAY(1..nn,datatype=integer[8]),
    XX::ARRAY(1..nn,datatype=integer[8]),
    shift::integer[4],
    pp::integer[8],
    LIB=lib
):

VandermondeSolve1 := proc( v::{Vector,list}, m::{Vector,list}, p::prime, shift::integer:=0 )
local t,i,j,M,x,a,L,R,y;

   t := numelems(v);
   printf("Quadratic C code Vandermonde solver: t=%d  p=%d\n",t,p);
   if numelems(m) <> t then error "v and m must be the same size"; fi;

   R := Array(0..t-1,datatype=integer[8]);
   y := Array(0..t-1,datatype=integer[8]);
   for i from 1 to t do R[i-1] := m[i]; y[i-1] := v[i]; od;

   a := Array(0..t-1,datatype=integer[8]);
   M := Array(0..t,datatype=integer[8]);

   VSolve1(R,y,t,a,M,shift,p);
   return [seq( a[i], i=0..t-1 )];

end:

t := 100:
p := 2^31-1;
alpha := numtheory[primroot](p);
r := rand(1000000):
e := {seq( r(), i=1..t )}:
while nops(e)<t do e := e union {r()} od:
m := [seq( alpha &^ e[i] mod p, i=1..t )]:
c := rand(p):
f := add( c()*x^e[i], i=1..t ):
v := [seq( Eval(f,x=modp(alpha &^ i,p)) mod p, i=0..t-1 )]:
a := VandermondeSolve1( v, m, p ):
save a, "aoutput";
add( a[i]*x^e[i], i=1..t ) - f;