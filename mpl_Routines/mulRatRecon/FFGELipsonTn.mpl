# kernelopts(numcpus=8);

FFGE := proc(A::Matrix, b::Vector, Y::list(name))
#
# Solve Ax=b for x using fraction-free Gaussian elimination and
# and Lispon's fraction-free back substitution
# Output x, y and det(A) where x[i] = y[i]/det(A)
#
local n,m,B,mu,i,j,k,det,x,y,num,r,numterms,f,g,h,mons;
   numterms := proc(f) if f=0 then 0 elif type(f,`+`) then nops(f) else 1 fi end;
   n,m := op(1,A);
   if n<>m then error "Matrix must be square" fi;
   m := op(1,b);
   if m<>n then error "Matrix and vector must have the same dimension" fi;
   B := <A|b>;
   mu := 1:
   det := 1;
   for k to n-1 do
       i := k;
       while i<=n and B[i,k]=0 do i := i+1; od;
       if i>n then return 0 fi;
       if i>k then # interchange row i and k
          for j from k to n+1 do B[i,j],B[k,j] := B[k,j],B[i,j] od;
          det := -det;
       fi;
       for i from k+1 to n do
           for j from k+1 to n+1 do
               num := expand(B[k,k]*B[i,j]-B[i,k]*B[k,j]);
               divide(num,mu,evaln(B[i,j]));
               if i=k+1 and j=k+1 then lprint(i,numterms(num),numterms(B[i,j])) fi;
           od;
           B[i,k] := 0;
       od;
       mu := B[k,k];
   od;
   det := det*B[n,n];
   printf("#det=%d\n",numterms(det));
   # Lipson's back substitution
   y := Vector(n);
   y[n] := B[n,n+1];
   printf("#y[%d]=%d\n",n,numterms(B[n,n+1]));
   for i from n-1 by -1 to 1 do
       num := expand(B[i,n+1]*B[n,n]-add(B[i,j]*y[j],j=i+1..n));
       divide(num,B[i,i],evaln(y[i]));
       printf("#N[%d]=%d  #y[%d]=%d\n",i,numterms(num),i,numterms(y[i]));
   od;
   f := Vector(n);
   g := Vector(n);
   for i from 1 to n do
       h := gcd(y[i],B[n,n],evaln(f[i]),evaln(g[i]));
       printf("#f[%d]=%d #g[%d]=%d #h=%d",i,numterms(f[i]),
            i,numterms(g[i]),numterms(h));
       printf("   deg(f[%d])=%d deg(g[%d])=%d\n",i,degree(f[i]),i,degree(g[i]));
   od;
   r := {seq(Y[i]=ithprime(i),i=1..nops(Y))};
   for i from 1 to n do
       coeffs(f[i],Y,'mons');
       mons := subs(r,[mons]);
       m := max(op(mons));
       printf("f[%d]: max m[%d] = %d = %a\n",i,i,m,ifactor(m));
       coeffs(g[i],Y,'mons');
       mons := subs(r,[mons]);
       m := max(op(mons));
       printf("g[%d]: max m[%d] = %d = %a\n",i,i,m,ifactor(m));
   od;
   return det,f,g;
end:

with(LinearAlgebra):

A := ToeplitzMatrix([y1,y2,y3,y4],symmetric);
b := <1,1,1,1>;
d,f,g := FFGE(A,b,[y1,y2,y3,y4]);
factor(d);

(* 
X := [x1,x2,x3,x4,x5,x6]:
A := ToeplitzMatrix(X,symmetric):
b := <1,1,1,1,1,1>:
FFGE(A,b,X):

profile(ToeplitzMatrix,gcd,divide,expand);
X := [x1,x2,x3,x4,x5,x6,x7]:
A := ToeplitzMatrix(X,symmetric):
b := <1,1,1,1,1,1,1>:
FFGE(A,b,X):

unprofile(ToeplitzMatrix,gcd,divide,expand);
profile(ToeplitzMatrix,gcd,divide,expand);

X := [x1,x2,x3,x4,x5,x6,x7,x8]:
A := ToeplitzMatrix(X,symmetric):
b := <1,1,1,1,1,1,1,1>:
FFGE(A,b,X):
showprofile();

unprofile(ToeplitzMatrix,gcd,divide,expand);
profile(ToeplitzMatrix,gcd,divide,expand);
printf("\nn=9\n");
X := [x1,x2,x3,x4,x5,x6,x7,x8,x9]:
A := ToeplitzMatrix(X,symmetric):
b := <1,1,1,1,1,1,1,1,1>:
FFGE(A,b,X):
showprofile();

unprofile(ToeplitzMatrix,gcd,divide,expand);
profile(ToeplitzMatrix,gcd,divide,expand);
printf("\nn=10\n");
X := [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10]:
A := ToeplitzMatrix(X,symmetric):
b := <1,1,1,1,1,1,1,1,1,1>:
FFGE(A,b,X):
showprofile();

unprofile(ToeplitzMatrix,gcd,divide,expand);
profile(ToeplitzMatrix,gcd,divide,expand);
printf("\nn=11\n");
X := [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11]:
A := ToeplitzMatrix(X,symmetric):
b := <1,1,1,1,1,1,1,1,1,1,1>:
FFGE(A,b,X):
showprofile();
quit;

unprofile(ToeplitzMatrix,gcd,divide,expand);
profile(ToeplitzMatrix,gcd,divide,expand);
printf("\nn=12\n");
X := [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]:
A := ToeplitzMatrix(X,symmetric):
b := <1,1,1,1,1,1,1,1,1,1,1,1>:
FFGE(A,b,X):
showprofile();

*)
