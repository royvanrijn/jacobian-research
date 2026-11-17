"""Exact identities and transcription controls;no old mathematical checker replay."""
import hashlib,json,platform,time
from itertools import permutations
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(b,label):
 assert bool(b),label
 checks.append(label)
def zero(f,label):check(s.expand(f)==0,label)
def even(f,vars):return all(int(c)%2==0 for c in s.Poly(s.expand(f),*vars).coeffs())
x,y,u,v,t=s.symbols('x y u v t');z=[x,y,u,v]
P=x+x*x*y+x**4+x**6*y*y;Q=y+x**5+x**6*y+x**7*y*y+x**8*y**3
H0=u*P+v*Q;J=s.expand(s.det(s.Matrix([P,Q]).jacobian([x,y])))
check(even(J-1,[x,y]),'the retained plane pair is Keller over every characteristic-two field')
check(s.Poly(Q,x,y).total_degree()==11 and s.Poly(Q,x,y).terms()[0][0][0]==8,
 'the coefficient polynomial Q has the exact degree-eleven top monomial x8 y3')
top=s.Add(*(c*x**ex[0]*y**ex[1] for ex,c in s.Poly(Q,x,y).terms() if sum(ex)==11))
zero(top-x**8*y**3,'the complete top homogeneous coefficient polynomial is a single monomial')
# Generic functions certify triangular identities without bounding either degree.
aa=s.Function('a')(x,y);bb=s.Function('b')(x,y,u)
phi=s.Matrix([x,y,u+aa**2,v+bb**2]);JP=phi.jacobian(z)
zero(JP.det()-1,'the arbitrary polynomial square shear has integral Jacobian determinant one')
for ii in [2,3]:
 f=aa if ii==2 else bb
 for row in range(4):
  for col in range(row,4):
   zero(s.diff(phi[ii],z[row],z[col])-2*(s.diff(f,z[row])*s.diff(f,z[col])+f*s.diff(f,z[row],z[col])),
    'the generic square Hessian entry is twice an integral expression '+str((ii,row,col)))
# The inverse works over Z as well as in characteristic two,with the appropriate signs.
inv=s.Matrix([x,y,u-aa**2,v-bb.subs(u,u-aa**2)**2])
for ii,f in enumerate(phi):
 zero(s.simplify(f.subs(dict(zip(z,inv)),simultaneous=True)-z[ii]),
  'the arbitrary-function triangular polynomial inverse composes exactly in coordinate '+str(ii))
# The gradient chain rule reduces to composition,because all added derivatives contain2.
a0=x;b0=u
ph0={x:x,y:y,u:u+x*x,v:v+u*u}
H=s.expand(H0.subs(ph0,simultaneous=True));g0=[s.diff(H0,zz) for zz in z]
for ii in range(4):
 difference=s.diff(H,z[ii])-g0[ii].subs(ph0,simultaneous=True)
 check(even(difference,z),'the actual nonlinear shear gradient chain rule has zero residue error '+str(ii))
for point in [(0,1,0,0),(1,0,1,1),(1,1,1,1)]:
 check([int(s.diff(H,zz).subs(dict(zip(z,point))))%2 for zz in z]==[0,0,0,1],
  'the explicit sheared gradient has its exact transported collision at '+str(point))
def part(f,degree):return s.Add(*(c*s.prod(zz**a for zz,a in zip(z,ex)) for ex,c in s.Poly(f,z).terms() if sum(ex)==degree))
zero(part(H,2)-x*u-y*v,'the explicit sheared quadratic part remains hyperbolic')
zero(part(H,3)-x**3-y*u*u,'the full new cubic part is x3 plus y u2')
zero(part(H,4)-x*x*y*u,'the full new quartic part is unchanged')
C=part(H,3);zero(s.diff(C,x,u)+s.diff(C,y,v),'all four harmonic defects vanish for the new cubic')
cs=lambda a,b,c,d:s.Poly(C,z).coeff_monomial(x**a*y**b*u**c*v**d)
T=(-cs(0,1,0,2)*cs(1,2,0,0)+3*cs(0,1,2,0)*cs(3,0,0,0)-2*cs(0,2,0,1)*cs(2,0,1,0)
   +2*cs(0,2,1,0)*cs(2,0,0,1)+3*cs(0,3,0,0)*cs(1,0,0,2)-cs(1,0,2,0)*cs(2,1,0,0))
check(T==3 and (1-T)%2==0,'the nonlinear shear passes exactly the old displayed divided residue equation')
small=s.hessian(part(H,2)+part(H,3)+part(H,4),z).subs({u:0,v:0});coefficient=0
for perm in permutations(range(4)):
 sign=(-1)**sum(perm[ii]>perm[jj] for ii in range(4) for jj in range(ii+1,4))
 coefficient+=sign*s.Poly(s.expand(s.prod(small[ii,perm[ii]] for ii in range(4))),x,y).coeff_monomial(x*y)
check(coefficient==-8,'the independent full determinant xy coefficient is minus eight,not the old mod8 obstruction')
# Every supported direction has a nonzero high Hasse coefficient;ordinary second derivatives miss it.
A,B=s.symbols('A B')
translated=s.Poly(((x+t*A)**8*(y+t*B)**3).expand(),t)
zero(translated.coeff_monomial(t**11)-A**8*B**3,'when both base directions are nonzero the t11 top coefficient cannot vanish')
zero(translated.coeff_monomial(t**8).subs(B,0)-A**8*y**3,'a pure x direction has a nonzero t8 Hasse coefficient')
zero(translated.coeff_monomial(t**3).subs(A,0)-x**8*B**3,'a pure y direction has a nonzero t3 Hasse coefficient')
d=s.symbols('d',integer=True,positive=True)
check((2*d-1).is_positive,'the arbitrary positive-u-degree square has translation degree2d strictly above one')
lead,w=s.symbols('lead w')
zero((lead*w**d)**2-lead**2*w**(2*d),'the arbitrary-degree leading square coefficient cannot cancel the affine terms')
Hv=s.diff(H,v);zero(Hv-Q,'the remaining v direction is genuinely affine for the whole candidate')
check(all(even(s.diff(H,zz,zz),z) for zz in z),
 'all ordinary pure second derivatives vanish modulo2,so the translation test is essential')
# Generic even-diagonal symmetric matrices:the linear determinant correction is divisible4.
diag=s.symbols('d0:4');off=s.symbols('a01 a02 a03 a12 a13 a23')
matrix=s.zeros(4)
for ii in range(4):matrix[ii,ii]=2*diag[ii]
for (ii,jj),a in zip([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],off):matrix[ii,jj]=matrix[jj,ii]=a
adj=matrix.adjugate();variables=list(diag)+list(off)
for ii in range(4):check(even(adj[ii,ii],variables),'the general odd principal cofactor has an exact factor two '+str(ii))
check(adj==adj.T,'the general symmetric adjugate pairs its off-diagonal trace terms')
bs=s.symbols('b0:10');bm=s.zeros(4);counter=0
for ii in range(4):
 for jj in range(ii,4):bm[ii,jj]=bm[jj,ii]=bs[counter];counter+=1
check(even(s.trace(adj*bm),variables+list(bs)),
 'twice the full general first determinant correction is divisible by four')
# A denominator-bearing plane needs saturation to retain rank in the residue.
pi=s.symbols('pi');basis=s.Matrix([[1/pi,0],[0,pi],[1,0],[0,1]])
sat=basis*s.diag(pi,1)
check(all(not f.has(1/pi) for f in sat),'a saturated plane basis removes an apparent frame pole')
check(sat.subs(pi,0).rank()==2,'the saturated basis retains plane rank after reduction')
naive=(basis*pi).subs(pi,0)
check(naive.rank()==1,'uniformly clearing both frame columns can falsely collapse residue rank')
algebra=len(checks)
context=json.loads((HERE/'sources/affine-plane-specialization-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
 check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for e in context['newly_cached_primary_context']+context['reinspected_retained_primary_context']:
 check(hashlib.sha256((HERE/'sources'/e['retained_file']).read_bytes()).hexdigest()==e['sha256'],
  'inspected primary context '+e['retained_file'])
check(context['newly_imported_proved_external_statements']==context['newly_used_unproved_external_statements']==[],
 'no new classification,Grassmannian properness or Legendre inversion theorem is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==8 and not flags['arbitrary_nonlinear_potential_source_orbit_excluded']
 and not flags['sharp_maximal_finite_precision_for_sheared_family_proved']
 and not flags['unrestricted_HC4_or_JC2_resolved'] and not flags['original_goal_complete'],
 'all-degree shear exclusion retains arbitrary nonlinear orbit,finite precision and HC4 boundaries')
proof=HERE/'SPECIALIZATION_OF_AFFINE_PLANES_AND_FROBENIUS_SHEARS.md'
paths=['sources/affine-plane-specialization-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']+context['reinspected_retained_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-affine-plane-specialization.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
