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
sys,xVar,yVar,nopsX,nopsY := genSystem("mike"):
sys;
xVar;
yVar;
nopsX;
nopsY;

(* Construct a black box B for the system of equations. *)
bkBox:= getBB(sys,xVar,yVar):

(* Choose a prime. *)
p := 107:
p;

solX := bkBox([1,2,3,4,5],p):
solX;

try
    num,den := MRFI(bkBox,nopsX,nopsY,yVar,p):
    catch:
    lprint("ERROR: ",lasterror()):
end try:

num;
den;
