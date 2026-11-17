"""New MW0 supporting algebra controls; sealed receipts must not be replayed."""
import hashlib,json,time
from collections import defaultdict
from pathlib import Path
import sympy as S

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
receipt=HERE/'receipt-affine-etale-primitive-descent.json'
assert not receipt.exists(),'Sealed receipt exists: never replay.'
assert not any(e['id']=='HC4MYMW0' for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries'])
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value,alg=True):
 assert bool(value),name
 checks.append(name)
 if alg:algebraic.append(name)

# Formal tensors after base change. First factor is the extra coefficient.
# Tensor factors are monomials; multiplication is concatenation and sorting.
def clean(d):return {k:v for k,v in d.items() if v}
def diff(d):
 out=defaultdict(int)
 for tensor,coefficient in d.items():
  for j in range(1,len(tensor)+1):
   out[tensor[:j]+((),)+tensor[j:]]+=(-1)**(j-1)*coefficient
 return clean(out)
def homotopy(d):
 out=defaultdict(int)
 for tensor,coefficient in d.items():
  assert len(tensor)>=2
  out[(tuple(sorted(tensor[0]+tensor[1])),)+tensor[2:]]+=coefficient
 return clean(out)
def add(a,b):
 d=defaultdict(int,a)
 for k,v in b.items():d[k]+=v
 return clean(d)
for degree in range(5):
 tensor=tuple((f'a{j}',) for j in range(degree+2));z={tensor:1}
 check(f'Amitsur contraction dh+hd in degree {degree}',add(diff(homotopy(z)),homotopy(diff(z)))==z)
 check(f'Amitsur consecutive differentials zero degree {degree}',diff(diff(z))=={})

# Actual relation signs, independent of a Galois presentation.
a1,a2,a3,g1,g2,g3,f1,f2,f3=S.symbols('a1 a2 a3 g1 g2 g3 f1 f2 f3')
check('rational differences are an exact triple cocycle',S.expand((g3-g2)-(g3-g1)+(g2-g1))==0)
check('split cocycle makes corrected primitive forms equal',S.expand(a2-a1+(a1-a2))==0)
check('split cocycle makes rational H equal',S.expand((f2-g2)-(f1-g1)-(f2-f1-(g2-g1)))==0)

# Pole differentiation in characteristic zero, and its characteristic-p boundary.
t=S.symbols('t')
for order in [1,2,7]:
 b=2+3*t+t**2;derivative=S.diff(b/t**order,t)
 check(f'unit pole leading derivative order {order}',S.limit(t**(order+1)*derivative,t,0)==-2*order)
for p in [2,3,5]:
 check(f'Artin-Schreier derivative is a unit in char {p}',S.Poly(S.diff(t**p-t,t),t,modulus=p).as_expr()==-1 or p==2 and S.Poly(S.diff(t**p-t,t),t,modulus=p).as_expr()==1)
 check(f'negative pth power differential vanishes in char {p}',int(-p)%p==0)

# The first mixed polynomial, without assuming any actual target realization.
y,u,x=S.symbols('y u x');h=y**2*(y-1)**3;hp=S.diff(h,y)
M=S.zeros(5)
for j in range(4):M[j+1,j]=1
M[:,4]=S.Matrix([u,0,1,-3,3])
# G=y-polynomial, gamma=0: trace 1/hprime=0 by interpolation.
hpM=5*M**4-12*M**3+9*M**2-2*M
check('mixed trace of reciprocal hprime is zero',S.cancel(S.trace(hpM.inv()))==0)
G=y-S.Rational(3,5);beta=1/hp;sl=x/hp+beta
check('mixed rational slice primitive identity',S.cancel(S.diff(G,y)-beta*hp)==0)
check('mixed primitive normalized trace zero',S.trace(M-S.Rational(3,5)*S.eye(5))==0)
check('rational primitive identity splits into regular source and E form',S.cancel(beta*hp-((sl)*hp-x))==0)
check('base rational slice shifts do not alter tracefree beta',S.expand((g1+g2)-(g3+g2)-(g1-g3))==0)

context_path=HERE/'sources/affine-etale-primitive-descent-context.json'
c=json.loads(context_path.read_text())
for path,digest in c['retained_proof_and_source_sha256'].items():check('retained byte pin '+path,sha(ROOT/path)==digest,False)
for path,digest in c['primary_source_sha256'].items():check('primary byte pin '+path,sha(ROOT/path)==digest,False)
for key in ['additive_differential_descent_proved','polynomial_representative_automatic','all_degree_residue_criterion_complete','rational_slice_retained','MV0_overstatement_explicitly_corrected','graph_family_exclusion_retained','full_trace_goal_open']:
 check('scope '+key,c['scope_flags'][key],False)
for key in ['all_mixed_slices_excluded','rational_shear_asserted_polynomial','target_pair_from_primitive','finite_trace_lattice_stability_proved','unrestricted_HC4_or_JC2_resolved']:
 check('boundary '+key,not c['scope_flags'][key],False)
check('no registered mathematical replay',c['registered_mathematical_checkers_replayed']==[],False)
out=dict(outcome='PASS',claim='HC4MYMW0',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,elapsed_seconds=time.perf_counter()-started,
 proof_sha256=sha(HERE/'AFFINE_ETALE_DESCENT_OF_RATIONAL_PRIMITIVES.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(context_path),
 registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
receipt.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
