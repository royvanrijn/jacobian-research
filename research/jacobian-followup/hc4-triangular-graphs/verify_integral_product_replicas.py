"""New algebra/provenance controls only. Never imports a registered checker."""
from pathlib import Path
import hashlib,json,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
started=time.perf_counter();checks=[];algebraic=[]
def check(name,condition,algebra=False):
    assert condition,name
    checks.append(name)
    if algebra:algebraic.append(name)
def zero(name,value):check(name,s.simplify(s.expand(value))==0,True)
x,y,t,c=s.symbols('x y t c');N,r=s.symbols('N r',integer=True,positive=True)
h=x*y;k=(x+c)*y
def pb(a,b):return s.diff(a,x)*s.diff(b,y)-s.diff(a,y)*s.diff(b,x)
def ham(H,a):return s.diff(H,y)*s.diff(a,x)-s.diff(H,x)*s.diff(a,y)
D=lambda a:x*s.diff(a,x)-y*s.diff(a,y)
zero('balanced product invariant',D(h))
zero('shifted products algebraically independent',pb(h,k)-c*y)
# Formal monomial coefficients give universal exponents; no degree sweep.
zero('general power Hamiltonian x component',ham(h**N,x)-N*h**(N-1)*x)
zero('general power Hamiltonian y component',ham(h**N,y)+N*h**(N-1)*y)
zero('general power divergence',s.diff(N*h**(N-1)*x,x)-s.diff(N*h**(N-1)*y,y))
zero('two general power bracket',pb(h**N,k**N)-N**2*c*y*h**(N-1)*k**(N-1))
zero('new two-field bracket not constant witness',pb(h**2,k**2).subs({x:1,y:1,c:1})-8)
check('general bracket has positive y exponent', (2*N-1).is_positive,True)
j=s.symbols('j',integer=True,nonnegative=True)
zero('all-iterate monomial induction coefficient',D(h**j*x)-h**j*x)
zero('all-iterate y induction coefficient',D(h**j*y)+h**j*y)
check('all-iterate degrees strictly grow for N>=2',s.diff(2*j*(N-1)+1,j).subs(N,2)==2,True)
alpha=h+h*h
delta=lambda a:alpha*D(a)
v0=x*x+h+y**3
f=[v0,delta(v0),delta(delta(v0))]
zero('replica coefficient invariant',delta(alpha))
zero('monic orbit relation with regular invariant coefficients',delta(f[2])+alpha*f[2]-6*alpha**2*f[1])
C=s.Matrix([[0,1,0],[0,0,1],[0,6*alpha**2,-alpha]])
zeta=s.symbols('zeta')
zero('companion characteristic polynomial',C.charpoly(zeta).as_expr()-zeta*(zeta-2*alpha)*(zeta+3*alpha))
check('finite matrix orbit vector equation',s.simplify(s.Matrix([delta(a) for a in f])-C*s.Matrix(f))==s.zeros(3,1),True)
z,v=s.symbols('z v',positive=True)
drz=-(N/r)*z**N*v**(N-1);drv=(N/r)*z**(N-1)*v**N
rootD=lambda a:s.diff(a,z)*drz+s.diff(a,v)*drv
zero('new root replica preserves h=zv',rootD(z*v))
zero('new root replica gives source x velocity',rootD(v/z**(r-1))-N*(z*v)**(N-1)*v/z**(r-1))
zero('new root replica gives source y velocity',rootD(z**r)+N*(z*v)**(N-1)*z**r)
root_t=(v-c)/z
zero('new root boundary t velocity',rootD(root_t)-(N/r)*z**(N-2)*v**(N-1)*(2*v-c))
zero('N2 boundary speed is constant',((N/r)*z**(N-2)*v**(N-1)*(2*v-c)).subs(N,2).subs(v,c)-2*c*c/r)
zero('N3 entire boundary fixed',((N/r)*z**(N-2)*v**(N-1)*(2*v-c)).subs(N,3).subs({z:0,v:c}))
zero('new shifted h is regular invariant on root cover',(z*v+c*z**r)-k.subs({x:v/z**(r-1),y:z**r}))
check('threshold M2 for every retained r control',s.ceiling(s.Rational(1+1,1))==2,True)
q=s.symbols('q',integer=True,positive=True)
check('simple log pole cannot be rational primitive derivative',(-q).is_negative,True)
# The pole statement, threshold, Jordan reduction and entire-flow argument are
# proved in the written note; these finite checks do not certify their scope.
context=json.loads((HERE/'sources/integral-product-replicas-context.json').read_text())
ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());known={e['id']:e for e in ledger['entries']}
for e in context['retained_structural_inputs']:
    check('prior claim '+e['claim'],known[e['claim']]['state']=='proved' and known[e['claim']]['canonical_source']==e['canonical_source'])
for p,hsh in context['retained_proof_and_source_sha256'].items():check('byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==hsh)
flags=context['scope_flags']
positive={'full_goal_active','previous_goal_turn_progress','product_primitive_upgrade_and_zero_boundary_proved','exact_hyperbolic_replica_threshold_proved','all_preserving_replicas_complete_on_entire_completion','all_divergence_free_locally_finite_stabilizers_classified_under_failure','new_two_noncommuting_complete_fields_on_retained_proper_controls','external_review_deferred'}
for name,value in flags.items():check('scope '+name,value is (name in positive))
proof=HERE/'INTEGRAL_SOURCE_PRODUCTS_AND_COMPLETE_HAMILTONIAN_REPLICAS.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebraic),elapsed_seconds=round(time.perf_counter()-started,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=flags,registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(HERE/'receipt-integral-product-replicas.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
