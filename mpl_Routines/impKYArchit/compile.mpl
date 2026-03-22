restart: 
with(LinearAlgebra):

read "getData.mpl":
read "mkBB.mpl":
read "helperFunc.mpl"
read "mqrfr.mpl":

(* Generate the system of equations. *)
sys,xVar,yVar,numY,numX := genSystem("bspline"):
sys;
xVar;
yVar;
numY;
numX;

(* Construct a black box B for the system of equations. *)
bkBox:= getBB(sys,xVar,yVar):

(* Choose a prime. *)
p := prevprime(2^63-1):

