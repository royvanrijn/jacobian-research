"""Finite algebraic premises for the written F1 adjoint obstruction."""
import hashlib,json,sys
from pathlib import Path
import sympy as S
t,s=S.symbols('t s');a=S.symbols('a0:4')
F=sum(a[i]*t**i for i in range(4))
q=[a[1]+a[2]*t+a[3]*t*t,a[2]+a[3]*t,a[3]]
r=[-a[0],-a[1]-a[0]*s,-a[2]-a[1]*s-a[0]*s*s]
for i in range(3):
 assert S.expand(t**i*(t*q[i]-r[i].subs(s,1/t))-F)==0
# At a finite common zero, q3=q2=q1=F=0 successively force every ai=0.
finite=S.Matrix([[S.expand(v).coeff(x) for x in a] for v in [F,*q]])
assert finite.det()==1
# At infinity the curve equation is a3=0; the other three rows are the
# regular transformed sections at s=0. Their common zero again forces ai=0.
infinity=S.Matrix([[S.expand(v).coeff(x) for x in a] for v in [a[3],*[v.subs(s,0) for v in r]]])
assert abs(infinity.det())==1
# F1 intersection form: (a C0+b F).(c C0+d F)=-ac+ad+bc.
def ip(v,w):return -v[0]*w[0]+v[0]*w[1]+v[1]*w[0]
R=(4,6);K=(-2,-3);L=(2,1);C0=(1,0)
assert 1+(ip(R,R)+ip(R,K))//2==9
assert ip(C0,R)==2 and ip(L,C0)==-1
h0=lambda aa,bb:sum(max(bb-j+1,0) for j in range(aa+1))
assert h0(2,1)==h0(1,1)==3
# On P1 x E, O(3) external-product L4 has square 24 and K-intersection -8.
assert 1+(2*3*4-2*4)//2==9
assert S.gcd(3,4)==1
assert 4*0+3*1+(4-1)*(3-1)==9
root=Path(__file__).resolve().parents[2]
result={'status':'PASS','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'finite_sections':[str(v) for v in q],'infinity_sections':[str(v) for v in r],'gluing_identities':'t^i*(t*q_i-r_i)=F for i=0,1,2','finite_coefficient_matrix_determinant':str(finite.det()),'infinity_coefficient_matrix_determinant':str(infinity.det()),'halving_genus':9,'F1_adjoint_h0':3,'fixed_divisor_degree':2,'product_arithmetic_genus':9,'geometric_argument_formally_verified':False,'cubic_points_enumerated':False,'correlated_target_complete':False}
packet=root/'artifacts/generated-results/elkies-k3-f1-trielliptic-obstruction-v1'
if '--write' in sys.argv:
 packet.mkdir(exist_ok=True);(packet/'result.json').write_text(json.dumps(result,indent=2)+'\n')
else:assert result==json.loads((packet/'result.json').read_text())
print('PASS: adjoint section gluing, basepoint coefficient matrices and intersection premises')
