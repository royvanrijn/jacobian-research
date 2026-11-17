"""New symbolic parameter transcriptions;never replay a registered checker."""
import hashlib,json,time
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];started=time.perf_counter();checks=[];algebra=[]
def check(name,test,alg=False):
    assert test,name
    checks.append(name)
    if alg:algebra.append(name)
def zero(name,z):check(name,s.simplify(s.expand(z))==0,True)
x,y,u,t,p,q=s.symbols('x y u t p q');n=s.symbols('n',integer=True,positive=True)
P=u*y-1;Q=y*P+u;du=n*u/x
zero('unit ideal',y*Q-(y*y+1)*P-1)
zero('fiber first coordinate elimination',(P-p).subs(u,q-p*y)+(p*y*y-q*y+p+1))
zero('fiber second coordinate implication',Q-q-y*(P-p)-(u+p*y-q))
elim=p*y*y-q*y+p+1
zero('generic quadratic discriminant',s.discriminant(elim,y)-(q*q-4*p*(p+1)))
zero('no identically zero fiber quadratic',elim.subs({p:0,q:0})-1)
zero('p0 fiber unique y',elim.subs(p,0).subs(y,1/q))
zero('p0 retained u nonzero',(q-p*y).subs(p,0)-q)
zero('zero u forces p minus one',((P-p).subs(u,0))+p+1)
Pu=s.diff(P,u);Py=s.diff(P,y);Qu=s.diff(Q,u);Qy=s.diff(Q,y)
zero('Jacobian in u y',Pu*Qy-Py*Qu-(y*P-u))
A=n*u/x;B=u*u-P*P
zero('pencil dx coefficient',Q*Pu*du-P*Qu*du-A)
zero('pencil dy coefficient',Q*Py-P*Qy-B)
zero('unimodular pencil coefficient identity',((2*x*y+x*u*(1-y*y))/n)*A-B-1)
zero('pencil derivative equals minus twice Jacobian',s.diff(B,u)*du+2*(n*u/x)*(y*P-u))
At=y-t*y*y-t;Bt=1-t*y
zero('pencil cyclic equation',P-t*Q-(u*At-Bt))
zero('numerator root denominator value',At.subs(y,1/t)+t)
check('numerator simple zero',s.diff(Bt,y)==-t,True)
zero('denominator discriminant',s.discriminant(At,y)-(1-4*t*t))
zero('infinite pencil cyclic equation',Q-(u*(y*y+1)-y))
check('infinite pencil simple numerator zero',s.diff(y,y)==1 and (y*y+1).subs(y,0)==1,True)
zero('t0 fiber Laurent identity',(P.subs(u,x**n)).subs(y,x**(-n)))
# The four divisor orders are parameter-independent,not a genus census.
orders=[1,1,-1,-1]
check('four cyclic divisor points total order zero',len(orders)==4 and sum(orders)==0,True)
check('every order fully ramified for any n',all(s.gcd(n,abs(v))==1 for v in orders),True)
zero('all n Riemann Hurwitz',-2*n+4*(n-1)-(2*(n-1)-2))
zero('affine first Betti number',2*(n-1)+3-1-2*n)
zero('affine Euler characteristic',2-2*(n-1)-3-(1-2*n))
# Known seed affine equivalence,with i^2=-1.
P1=P.subs(u,x);Q1=Q.subs(u,x);aa=s.I*y;bb=s.I*x
zero('Jelonek known first seed coordinate',aa*(1+aa*bb)-bb+s.I*Q1)
zero('Jelonek known second seed coordinate',aa*bb+P1+1)
# Universal coframe adjugate algebra at one formal jet.
px,py,qx,qy,pp,qq=s.symbols('px py qx qy pp qq');J=px*qy-py*qx
Ex=(pp*qy-qq*py)/J;Ey=(-pp*qx+qq*px)/J
zero('actual Euler P action',px*Ex+py*Ey-pp)
zero('actual Euler Q action',qx*Ex+qy*Ey-qq)
zero('actual contraction dx',-Ey+(qq*px-pp*qx)/J)
zero('actual contraction dy',Ex+(qq*py-pp*qy)/J)
zero('actual coframe inverse determinant',s.det(s.Matrix([[qy,-py],[-qx,px]]))/J**2-1/J)
ctx=json.loads((H/'sources/smooth-projective-pencils-context.json').read_text());known={e['id']:e for e in json.loads((R/'MATH_STATUS.json').read_text())['entries']}
check('new checker is not yet registered','HC4MYMD0' not in known)
for e in ctx['retained_structural_inputs']:check('prior '+e['claim'],known[e['claim']]['state']=='proved' and known[e['claim']]['canonical_source']==e['canonical_source'])
for path,d in ctx['retained_proof_and_source_sha256'].items():check('byte '+path,hashlib.sha256((R/path).read_bytes()).hexdigest()==d)
yes={'full_goal_active','previous_genus_controls_and_published_seed_credited','all_n_whole_plane_quasifinite_punctured_image_control','smooth_surjective_complete_base_projective_pencil','every_pencil_fiber_irreducible','generic_genus_n_minus_one_and_three_punctures','Jacobian_nonconstant_and_not_Keller','actual_missing_value_radial_field_has_divergence_two','constant_differential_mark_not_supplied_by_controls','Chau_regular_value_nonempty_convention_preserved','external_review_deferred'}
for path,value in ctx['scope_flags'].items():check('scope '+path,value is (path in yes))
proof=H/'SMOOTH_PROJECTIVE_PENCILS_WITH_ARBITRARY_GENUS.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-started,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=ctx['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-smooth-projective-pencils.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({p:out[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
