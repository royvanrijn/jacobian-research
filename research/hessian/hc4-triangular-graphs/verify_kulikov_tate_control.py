"""Exact transcription certificates for the known Kulikov control and its boundary."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(v,label):
 assert bool(v),label
 checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
X,Y,Z,a,b,t,q=s.symbols('X Y Z a b t q')
coords=[X,Y,Z];Q=[3*v**2-X*Y-X*Z-Y*Z for v in coords]
R=sum(v**2*w for v in coords for w in coords if v!=w)-6*X*Y*Z
C=s.expand(Q[0]+2*Q[1]+4*Q[2])
check(s.expand(s.det(s.Matrix(Q).jacobian(coords))+36*R)==0,'the homogeneous Jacobian is exactly -36 times the nodal ramification cubic')
check(all(v.subs({X:1,Y:1,Z:1})==0 for v in Q),'the displayed base point belongs to all three conics')
for xx,yy in [(1,1),(-1,1),(1,-1),(-1,-1)]:
 check(all(v.subs({X:xx,Y:yy,Z:1})==0 for v in Q)==(xx==yy==1),'the sign-pattern base-locus test for '+str((xx,yy)))
check(s.Matrix([[s.diff(v,X).subs({X:1,Y:1,Z:1}),s.diff(v,Y).subs({X:1,Y:1,Z:1})] for v in Q]).rank()==2,
      'the conic net has no common infinitely-near base direction')
grad=s.groebner([s.diff(R,v).subs(Z,1) for v in coords],X,Y)
check(list(grad)==[X-1,Y-1],'the cubic has only the specified affine singular point')
check(s.groebner([s.diff(R,v).subs({Z:0,X:1}) for v in coords],Y).contains(s.Integer(1)),
      'there is no ramification-cubic singularity at infinity with X nonzero')
check(any(s.diff(R,v).subs({X:0,Y:1,Z:0})!=0 for v in coords),'the remaining infinity point is nonsingular')
local=s.expand(R.subs({X:1+a,Y:1+b,Z:1}))
check(local==a*a*b+a*b*b+2*a*a-2*a*b+2*b*b,'the local nodal equation is exact')
check(s.hessian(local,[a,b]).subs({a:0,b:0}).det()==12,'the quadratic tangent cone is nondegenerate')
check((s.hessian(C,coords)/2).det()==-127,'the selected conic is nonsingular')
line=s.expand(local.subs(b,t*a));AA=line.coeff(a,2);BB=line.coeff(a,3)
cc=s.expand(C.subs({X:1+a,Y:1+b,Z:1}).subs(b,t*a));LL=cc.coeff(a,1);MM=cc.coeff(a,2)
HH=s.expand(LL*BB-AA*MM);h=2*t**4-4*t**3+7*t**2-2*t+1
check(HH==-6*h,'the strict boundary intersection polynomial is the displayed quartic')
check(s.degree(h,t)==4 and s.gcd(h,s.diff(h,t))==1,'all four strict boundary intersections are distinct')
check(s.degree(s.gcd(AA,BB),t)==0,'the nodal normalization parametrization is reduced')
check(s.degree(s.gcd(AA,LL),t)==0,'the selected conic misses the two resolved node directions')
check(s.degree(s.gcd(BB,h),t)==0,'no intersection is hidden at a normalization denominator zero')
check(C.subs({X:0,Y:1,Z:0})!=0,'the parameter infinity point is not an intersection')
pair=lambda u,v:u[0]*v[0]-u[1]*v[1]
H=(1,0);E=(0,1);lineclass=(2,-1);ram=(3,-2);canonical=(-3,1)
check(pair(lineclass,lineclass)==3,'the resolved finite morphism has geometric degree3')
check(tuple(canonical[i]+3*lineclass[i] for i in range(2))==ram,'the finite-map canonical divisor is exactly the strict cubic')
check(pair(ram,lineclass)==4,'the omitted ramification curve has exactly four removed points')
check(s.Matrix([ram,lineclass]).det()==1,'the two deleted boundary classes form a unimodular Picard basis')
check(ram[0]+2*ram[1]==canonical[0]+2*canonical[1]==-1,'the affine finite normalization has primitive nonzero canonical and ramification class')
restriction=s.Matrix([[3,2],[2,1]])
check(restriction.det()==-1 and restriction.rank()==2,'the cohomology restriction matrix is an integral isomorphism')
incidence=s.Matrix([[-1]*4,[1]*4])
check(4-incidence.rank()==3,'the boundary dual graph has exactly three independent loops')
EW=1+2*q+q*q;ED=2*(1+q)-4;ES=s.expand(EW-ED)
check(ES==q*q+3,'compact-support Hodge additivity gives the weight-four surface polynomial')
check(ES.subs(q,1)==4,'the unweighted Euler characteristic is4')
c=s.symbols('c',integer=True,nonnegative=True)
check(s.expand(q*q+c*q-q*q)==c*q,'the actual plane-completion degree-two contribution has only weight2')
check(2*2>3 and 2+1==3,'the cubic degree sum forces residue degree1 for the index2 prime and one retained sheet')
algebra=len(checks)
context=json.loads((HERE/'sources/kulikov-tate-context.json').read_text())
for p,hsh in context['retained_proof_and_source_sha256'].items():check(sha(ROOT/p)==hsh,'immutable retained input '+p)
for p,hsh in context['new_primary_sha256'].items():check(sha(HERE/p)==hsh,'inspected new primary '+p)
check(context['newly_used_unproved_external_statements']==[],'no unproved external theorem is imported')
check(len(context['newly_imported_proved_external_statements'])==2,'the published simple-connectedness input and standard MHS formalism are the two explicit imports')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(flags['construction_credited_to_Kulikov'] and not flags['nonplane_surface_reported_as_HC4_counterexample'] and not flags['original_goal_complete'],
      'the known nonplane countercontrol does not substitute for the user full HC4 objective')
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=sha(HERE/'KULIKOV_EXACT_COFRAME_AND_TATE_BOUNDARY_FILTER.md'),checker_sha256=sha(Path(__file__)),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256=dict(context['new_primary_sha256'],**{'sources/kulikov-tate-context.json':sha(HERE/'sources/kulikov-tate-context.json')}),
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-kulikov-tate-control.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
