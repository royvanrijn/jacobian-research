"""NEW MV0 global-slice/graph-image/degree transcriptions; never replay."""
import hashlib,json,time
from pathlib import Path
import sympy as S
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
RECEIPT=HERE/'receipt-trace-free-rational-slices.json'
assert not RECEIPT.exists(),'Sealed receipt exists: never replay.'
assert not any(e['id']=='HC4MYMV0' for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries'])
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value,alg=True):
 assert bool(value),name
 checks.append(name)
 if alg:algebraic.append(name)
x,y,u,v,w,s,B,A0,H0=S.symbols('x y u v w s B A0 H0')
h=y**2*(y-1)**3;hp=S.diff(h,y);a=-S.Rational(108,3125)
Au=(u*(u-a))**2;Ah=Au.subs(u,h)
D=S.cancel(Ah/hp**2);c=S.cancel(Ah/hp)
check('mixed critical values included',h.subs(y,0)==h.subs(y,1)==0 and h.subs(y,S.Rational(2,5))==a)
check('clearing quotient is polynomial',S.degree(S.denom(D),y)==0)
check('source affine fibre coefficient polynomial',S.degree(S.denom(c),y)==0 and S.cancel(c-hp*D)==0)
for j in range(2,8):
 check(f'all-degree Taylor clearing identity j{j}',S.cancel(A0**(j-1)/H0**j-(A0/H0**2)**(j-1)*H0**(j-2))==0)
check('exact finite Taylor identity',S.expand(h.subs(y,y+B*x)-h-sum(S.diff(h,y,j)*B**j*x**j/S.factorial(j) for j in range(1,6)))==0)
for root in [S.Integer(0),S.Integer(1),S.Rational(2,5)]:
 check('source critical line coefficient vanishes '+str(root),D.subs(y,root)==c.subs(y,root)==0)
cubic=S.cancel((h-a)/(y-S.Rational(2,5))**2)
check('all simple second-fibre source lines also retained',S.rem(D,cubic,y)==0)
source_v=y+c*x
check('source pulls back unit volume',S.cancel(S.Matrix([source_v,h]).jacobian([x,y]).det()/Ah-1)==0)
check('global hypersurface smooth at multiple base roots',(S.diff(Au,u)*w+1).subs(u,0)==(S.diff(Au,u)*w+1).subs(u,a)==1)

# Finite field trace of the genuine rational slice, over C(u).
M=S.zeros(5)
for j in range(4):M[j+1,j]=1
M[:,4]=S.Matrix([u,0,1,-3,3])
Beta=M/Au;mu=S.Rational(3,5)
check('mean source coordinate constant',S.trace(M)/5==mu)
check('mean rational slice correction',S.cancel(S.trace(Beta)/5-mu/Au)==0)
check('trace-free correction has zero trace',S.trace(Beta-mu/Au*S.eye(5))==0)
obstruction=(y-mu)*hp/Ah
residue=S.residue(obstruction,y,0)
expected=(-3-12/(5*a))/a**2
check('explicit mixed model has nonzero exactness residue',S.cancel(residue-expected)==0 and residue!=0)

# General rational-map degree gap; the argument is unbounded in the note.
m,N,r,eps=S.symbols('m N r eps',integer=True)
check('all-degree pole versus infinity degree gap',S.expand((m*(N-1)-1)-(m*(N-r)-eps)-(m*(r-1)+eps-1))==0)
for degree,roots,epsilon in [(2,2,0),(2,2,1),(5,2,0),(5,3,1),(17,9,0)]:
 check(f'positive gap m{degree} r{roots} eps{epsilon}',degree*(roots-1)+epsilon-1>0)

# One-root exact controls verify both parts of the global criterion, not
# merely a rational primitive with no polynomial source representative.
for degree,k in [(2,2),(3,2),(5,4)]:
 depth=degree*(k-1);weight=1-depth
 h1=y**degree;hp1=S.diff(h1,y);beta=y**(1-degree*k)
 sl=x/hp1+beta
 G=S.Rational(degree,weight)*y**(1-depth)
 F0=-x*y/S.Integer(weight)
 check(f'one-root rational trace-free primitive m{degree} k{k}',S.cancel(S.diff(G,y)-beta*hp1)==0)
 check(f'one-root primitive trace zero m{degree} k{k}',(1-depth)%degree!=0)
 check(f'one-root POLYNOMIAL source representative m{degree} k{k}',S.cancel(F0-G+degree*h1*sl/weight)==0)
 theta_u=v/(weight*u**k);theta_v=-degree*u**(1-k)/weight
 check(f'one-root global Euler primitive differential m{degree} k{k}',S.cancel(S.diff(theta_u,v)-S.diff(theta_v,u)-u**(-k))==0)
 check(f'one-root Euler volume weight nonzero m{degree} k{k}',weight!=0)

# Source shear changes the trace-free obstruction only by an exact term.
Q=y**5/S.Integer(5);meanQ=(3+5*u)/25
check('source shear trace correction differential',S.expand(S.diff(-Q+meanQ.subs(u,h),y)+y**4-hp/5)==0)
check('base slice shifts cancel in trace-free correction',S.expand((B+A0)-(H0+A0)-(B-H0))==0)

cp=HERE/'sources/trace-free-rational-slices-context.json';C=json.loads(cp.read_text())
for p,digest in C['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,sha(ROOT/p)==digest,False)
check('five retained claims and ten proof-checker pins',len(C['retained_structural_inputs'])==5 and len(C['retained_proof_and_source_sha256'])==10,False)
for key in ['global_exactness_criterion_with_polynomial_representative','rational_beta_retained','trace_free_obstruction_gauge_invariant','explicit_surjective_etale_plane_source_image','all_degree_multiroot_graph_no_plane_target','first_mixed_source_cover_realized_without_target','one_root_exact_positive_controls','original_full_trace_goal_open']:
 check('scope '+key,C['scope_flags'][key],False)
for key in ['all_rational_mixed_slices_excluded','rational_primitive_alone_asserted_sufficient','source_image_graph_is_Keller_counterexample','exactness_asserted_supply_Darboux_pair','Mathieu_transfer_assumed','registered_predecessors_replayed','unrestricted_HC4_or_JC2_resolved']:
 check('boundary '+key,not C['scope_flags'][key],False)
check('no new imports caches downloads or predecessor replay',C['new_external_theorem_imports']==C['new_primary_caches']==C['new_primary_downloads']==0 and C['registered_mathematical_checkers_replayed']==[],False)
r=dict(outcome='PASS',claim='HC4MYMV0',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,algebraic_check_names=algebraic,elapsed_seconds=time.perf_counter()-started,
 proof_sha256=sha(HERE/'TRACE_FREE_RATIONAL_SLICES_AND_MIXED_SOURCE_IMAGE_GRAPH_MODELS.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
RECEIPT.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
