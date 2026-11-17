needsPackage "BernsteinSato";
A=QQ[u,v]/ideal(u^2+v^2);K=frac A;
print("FIELD");print K;print(u^2+v^2);
W=K[x,y,t,dx,dy,dt,WeylAlgebra=>{x=>dx,y=>dy,t=>dt}];
I=ideal(dx+x-t*u,dy+y-t*v,dt-u*x-v*y);
M=Dintegration(0,I,{1,1,0});
print("SOURCE");print presentation M;
print("COMPLETE");exit 0;
