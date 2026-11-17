import sympy as S
from functools import lru_cache
from pathlib import Path
U,V,K,r=S.symbols('U V K r');a,b,c,d,e,A,B,C,D,E,lu,q0=S.symbols('a b c d e A B C D E lu q0')
m=-2*c*e**2/K**4
lv=2*D*c-2*b*d*m-d/(2*c)-d**3*e**2*m/(2*c**3)
# K remains an independent coordinate in differentiation; use root symbol kap.
kap=S.Symbol('kap');x=(U,V,kap)
F=[kap*U*V+e*U+c*V**2-m*kap**2/2,d*V**3-d*e*m*kap/(2*c),b*U**2-E*V*kap-q0,-D*kap*U-lv*V,a*U**3-B*V**2/2-lu*U,-C*U*V,-A*U**2/2]
p0={U:2*c*e/K**2,V:-e/K,kap:K}
g=[S.Matrix([S.diff(f,z).subs(p0) for z in x]) for f in F]
H=[S.hessian(f,x).subs(p0) for f in F]
def simp(z):return S.factor(z)
def smat(M):return M.applyfunc(simp)
Hi=smat(H[0].inv())
T=[[[[S.diff(f,x[i],x[j],x[k]) for k in range(3)] for j in range(3)] for i in range(3)] for f in F]
def bil(Tj,v,w):return S.Matrix([sum(Tj[i][j][k]*v[j]*w[k] for j in range(3) for k in range(3)) for i in range(3)])
def lin(Tj,v):return S.Matrix(3,3,lambda i,j:sum(Tj[i][j][k]*v[k] for k in range(3)))
z1=smat(-Hi*g[1]);print('z1',list(z1),flush=True)
z2=smat(-Hi*(bil(T[0],z1,z1)/2+H[1]*z1+g[2]));print('z2',list(z2),flush=True)
z3=smat(-Hi*(bil(T[0],z1,z2)+H[1]*z2+bil(T[1],z1,z1)/2+H[2]*z1+g[3]));print('z3 done',flush=True)
shift={v:p0[v]+r*z1[i]+r*r*z2[i] for i,v in enumerate(x)}
# Only monomials up to r^5 can contribute.
action5=0
for j in range(6):
 action5+=S.expand(F[j].subs(shift,simultaneous=True)).coeff(r,5-j)
action5=simp(action5);print('action5',action5,flush=True)
H1=smat(H[1]+lin(T[0],z1));H2=smat(H[2]+lin(T[1],z1)+lin(T[0],z2));H3=smat(H[3]+lin(T[1],z2)+lin(T[0],z3))
M1=smat(Hi*H1);M2=smat(Hi*H2);M3=smat(Hi*H3)
logdet3=simp(S.trace(M3)-S.trace(M1*M2)+S.trace(M1**3)/3);print('logdet3',logdet3,flush=True)
# Wick with covariance -H^-1; differentiate it by one dual variable s.
s=S.Symbol('s');cov=-Hi+s*(Hi*H1*Hi)
@lru_cache(None)
def wick(indices):
 if not indices:return S.Integer(1)
 first=indices[0];out=0
 for j in range(1,len(indices)):
  remaining=indices[1:j]+indices[j+1:]
  out+=cov[first,indices[j]]*wick(remaining)
 return S.expand(out)
L1=simp(S.expand(wick((0,0,1,1,2,2))).coeff(s,1)/2+d*wick((0,1,1,1,1,2)).subs(s,0))
print('stationary_L1',L1,flush=True)
chi3=simp(action5-logdet3/2+L1)
print('CHI3',chi3,flush=True)
Path('/home/royvanrijn/src/jacobian-research/research/.scratch/double_resonance_chi3.txt').write_text(str(chi3)+'\n')
