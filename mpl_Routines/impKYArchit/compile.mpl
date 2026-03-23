restart: 
with(LinearAlgebra):

read "getData.mpl":
read "mkBB.mpl":
read "helperFunc.mpl":
read "mqrfr.mpl":
read "ndsa.mpl":
read "getPtAffLn.mpl":

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
p := prevprime(2^63-1):
p;

solX := bkBox([1,2,3,4,5],p):
solX;

MQRFRres,linSys := NDSA(BB,[seq(1,i=1..2)],2,p,4,2):
MQRFRres;
linSys;
