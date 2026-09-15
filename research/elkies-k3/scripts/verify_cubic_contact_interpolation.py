"""Symbolic contact identities and controls, not an arithmetic K3 point search."""
import hashlib,json,sys
from pathlib import Path
import sympy as S
t=S.symbols('t'); c=S.symbols('c0:6'); z=S.symbols('z0:6')
mod=t**6+sum(c[i]*t**i for i in range(6)); U=sum(z[i]*t**i for i in range(6))
M=S.Matrix([[z[k],z[k-1]-z[5]*c[k]] for k in [3,4,5]])
assert all(S.expand(S.rem(t*U,mod,t)).coeff(t,k)==M[k-3,1] for k in [3,4,5])
minors=[S.expand(M[list(pair),:].det()) for pair in [(0,1),(0,2),(1,2)]]
f=t**3-t-1
assert S.Poly(f,t,modulus=2).is_irreducible

def reduce(v):return S.rem(v,f,t)
def inverse(v):return S.invert(v,f,t)
def check(u,v,expected_rank):
 u=reduce(u);v=reduce(v)
 lift=S.expand(u+f*reduce((v-S.diff(u,t))*inverse(S.diff(f,t))))
 assert reduce(lift-u)==0 and reduce(S.diff(lift,t)-v)==0
 columns=[lift,S.rem(t*lift,f*f,t)]
 matrix=S.Matrix([[S.expand(col).coeff(t,k) for col in columns] for k in [3,4,5]])
 assert matrix.rank()==expected_rank
 A,B,C,D,E=S.symbols('A B C D E'); unknowns=[A,B,C,D,E]
 value=reduce(A+B*t+C*t*t+(D+E*t)*u)
 slope=reduce(B+2*C*t+E*u+(D+E*t)*v)
 equations=[S.expand(q).coeff(t,k) for q in [value,slope] for k in range(3)]
 full=S.Matrix([[q.coeff(a) for a in unknowns] for q in equations])
 assert full.rank()==3+expected_rank
 out={'u':str(u),'v':str(v),'Hermite_lift':str(lift),'matrix':[[str(x) for x in row] for row in matrix.tolist()],'matrix_rank':expected_rank,'full_jet_rank':full.rank()}
 if expected_rank==1:
  den=matrix.nullspace()[0];den=S.expand(den[0]+den[1]*t)
  num=-S.rem(den*lift,f*f,t)
  assert S.degree(num,t)<=2 and S.gcd(den,f)==1 and S.gcd(num,den)==1
  # Homogeneous coprimality also requires no shared factor at infinity.
  assert S.degree(num,t)==2 or S.degree(den,t)==1
  assert S.rem(num+den*lift,f*f,t)==0
  out.update(numerator=str(-num),denominator=str(den))
 return out
u=reduce((t*t+t+1)*inverse(t+2));v=reduce((t*t+4*t+1)*inverse((t+2)**2))
controls=[check(u,v,1),check(u,v+1,2),check(t+1,1,0),check(t*t+t+3,2*t+1,1)]
assert S.cancel(S.sympify(controls[0]['numerator'])/S.sympify(controls[0]['denominator'])-(t*t+t+1)/(t+2))==0
# Intersection, genus, primitivity and the exact Castelnuovo--Severi boundary.
assert [(n,k) for n in range(3) for k in range(5) if -n+2*k==3]==[(1,2)]
assert -4+4*2+6==10
assert 10-2*3==4 and 4*(4+3)==28
assert 4*0+3*0+(4-1)*(3-1)==6<9
assert 4*0+3*1+(4-1)*(3-1)==9
assert max(3*(3-1)//2+1,3*1*(1-1)+1*(3*3-1-6))==4<9
root=Path(__file__).resolve().parents[2]
result={'status':'PASS','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'contact_matrix_minors':[str(m) for m in minors],'controls':controls,'arithmetic_genus':4,'quotient':'F1','section_self_intersection':3,'branch_intersection':10,'contact_degree':3,'residual_branch_degree':4,'anti_trace_height':28,'CS_genus0_bound':6,'CS_genus1_bound':9,'actual_halving_point_constructed':False,'trielliptic_map_constructed':False,'correlated_target_complete':False}
packet=root/'artifacts/generated-results/elkies-k3-cubic-contact-interpolation-v1'
if '--write' in sys.argv:
 packet.mkdir(exist_ok=True);(packet/'result.json').write_text(json.dumps(result,indent=2)+'\n')
else:assert result==json.loads((packet/'result.json').read_text())
print('PASS: cubic contact Hermite matrix, controls and genus premises; no K3 point or correlated pair')
