needsPackage "BernsteinSato";
W=QQ[x,dx,WeylAlgebra=>{x=>dx}];
I=ideal(dx);F=matrix{{x*dx}};
Q=F//gens I;
print("LIFT");print Q;print(gens I*Q);print(F);
print("PRODUCT");print(matrix{{dx}}*matrix{{x}});print(dx*x);
print("DONE");exit 0;
