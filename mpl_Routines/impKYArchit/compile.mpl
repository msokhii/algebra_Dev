restart: 
with(LinearAlgebra):

read "getData.mpl":
read "mkBB.mpl":
read "helperFunc.mpl":
read "mqrfr.mpl":
read "ndsa.mpl":
read "getPtAffLn.mpl":
read "mrfi.mpl":

(* Generate the system of equations. *)
sys,xVar,yVar,nopsX,nopsY := genSystem("bspline"):
sys;
xVar;
yVar;
nopsX;
nopsY;
numLines:=0:

(* Construct a black box B for the system of equations. *)
bkBox:= getBB(sys,xVar,yVar):

(* Choose a prime. *)
p := prevprime(2^63-1):
p;

solX := bkBox([1,2,3,4,5],p):
solX;

try
    num,den := MRFI(bkBox,nopsY,nopsX,yVar,p):
    catch:
    lprint("ERROR: ",lasterror()):
end try:

num;
den;
