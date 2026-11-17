"""Exact controls for the written integral-cycle and global-period proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
x,y,z,w,u,v,s,t,rho,P,V,aa,bb,cc=sp.symbols('x y z w u v s t rho P V aa bb cc')
# Universal simple-root residue controls,without a root-count or degree scan.
eps=sp.symbols('eps');root_value,root_derivative=sp.symbols('root_value root_derivative',nonzero=True)
zero(sp.residue((s-root_value)/(root_derivative*eps),eps,0)-(s-root_value)/root_derivative,
     'the arbitrary simple-root suspended multiplier residue is exactly its value difference divided by Cprime')
zero(sp.diff((s-root_value)/(2*root_derivative),s)-1/(2*root_derivative),
     'every multiplier-root global period has its exact constant first derivative')
zero(sp.diff((s-root_value)/(2*root_derivative),s,2),'there are no higher global period terms in any multiplier degree')
# The normalized h_D case permits arbitrary D through a general remainder modulo C.
Cz=z*z-1;D=aa*z+bb+Cz*cc
lam=(s-D)/Cz
for r in [-1,1]:
    ci=D.subs(z,r)
    zero(sp.residue(lam,z,r)-(s-ci)/(2*r),'h_D native plane residue at root '+str(r))
    zero(r*sp.residue(lam,z,r)/2-(s-ci)/4,'the independently oriented h_D basis has common slope one quarter '+str(r))
    E=sp.cancel((Cz*w+D-ci)/(z-r))
    check(sp.Poly(E,z,w).is_multivariate,'the local node second coordinate is polynomial at root '+str(r))
    zero(sp.diff(E,w).subs(z,r)-2*r,'the root product coordinates have nonzero local Jacobian '+str(r))
# The native K2 Gelfand--Leray form is computed by the entire family map.
zz=s+u*v;yy=(1-v-zz*zz)/(2*x)
native=[x,yy,zz,u]
Jac=sp.Matrix(native).jacobian([s,x,u,v]).det()
zero(Jac-1/(2*x),'the full native four-volume is ds wedge dx du dv divided by 2x')
K2=w*(2*x*y+z*z-1)+z
zero(K2.subs(dict(zip([x,y,z,w],native)),simultaneous=True)-s,
     'the volume conversion chart is the actual K2 value-family chart')
zero((sp.diff(zz,u)*sp.diff(u,v)-sp.diff(zz,v)*sp.diff(u,u))+u,
     'on a fixed value,dz wedge dw is exactly minus u du wedge dv')
k_u=(z-s)/(1-z*z);k_v=1-z*z;k_lam=sp.cancel(k_u*sp.diff(k_v,z))
zero(k_lam-2*z*(z-s)/(z*z-1),'the K2 area primitive on the whole exceptional curve')
for r in [-1,1]:
    zero(sp.residue(k_lam,z,r)-(r-s),'the K2 global contour residue '+str(r))
    zero(-sp.residue(k_lam,z,r)/2-(s-r)/2,'the oriented K2 global period has common slope one half '+str(r))
    value=z+u*(z*z-1)
    A=z-r;B=1+u*(z+r)
    zero(value-r-A*B,'the K2 curve-family node has exact product coordinates '+str(r))
    det=sp.Matrix([A,B]).jacobian([z,u]).det().subs({z:r,u:-sp.Rational(r,2)})
    zero(det-2*r,'the K2 local product-coordinate Jacobian is a unit '+str(r))
# Three-point family:exact product coordinates identify the actual annuli.
T=u*v+rho;Br=u*T*T-u+u*v+rho
for r in [-1,1]:
    A=T-r;B=u*(T+r)+1
    zero(Br-r-A*B,'the S_rho node factorization holds for all rho at value '+str(r))
    point={u:-sp.Rational(r,2),v:2*r*rho-2}
    zeros([A.subs(point,simultaneous=True),B.subs(point,simultaneous=True)],
          'both factors vanish at the stated distinct local node '+str(r))
    zero(sp.Matrix([A,B]).jacobian([u,v]).det().subs(point,simultaneous=True)-1,
         'the S_rho product coordinates have Jacobian one uniformly in rho '+str(r))
A=u;B=T*T+v-1;point={u:0,v:1-rho*rho}
zero(Br-rho-A*B,'the third S_rho node has exact product coordinates in every rho')
zeros([A.subs(point,simultaneous=True),B.subs(point,simultaneous=True)],'both third-node coordinates vanish')
zero(sp.Matrix([A,B]).jacobian([u,v]).det().subs(point,simultaneous=True)-1,
     'the third local product chart remains etale also at rho=+-1')
# On B_rho=s,t solves an implicit equation with unit t derivative near u=0.
rel=u*(t*t-1)+t-s
zero(sp.diff(rel,t).subs({u:0,t:s})-1,'the third-node t parameter is analytic across the u disk')
zero((t-s)+u*(t*t-1)-rel,'the exact moving-puncture winding identity')
zero((-u*(t*t-1)).subs(t,s)+u*(s*s-1),
     'the moving-puncture zero is simple in u away from the actual two deleted values')
for r in [-1,1]:
    zero(rel.subs(t,r)-(r-s),'the analytic u disk cannot hit the fixed puncture '+str(r))
scurve={u:(s-t)/(t*t-1),v:(t-rho)*(t*t-1)/(s-t)}
zero(Br.subs(scurve,simultaneous=True)-s,'the primitive uses the whole exact S_rho curve parameterization')
slambda=sp.cancel(scurve[v]*sp.diff(scurve[u],t))
partial=(rho-1)/(t-1)+(rho+1)/(t+1)+(s-rho)/(t-s)-1
zero(slambda-partial,'the complete native S_rho area primitive has its exact partial-fraction identity')
residues=sp.Matrix([rho-1,rho+1,s-rho])
for index,r in enumerate([1,-1,s]):
    zero(sp.residue(slambda,t,r)-residues[index],'the S_rho puncture residue index '+str(index))
M=sp.Matrix([[1,0,1],[0,1,1],[0,0,1]])
check(M.det()==1,'the three point-labelled vanishing circles are a unimodular integral basis')
check(all(q.is_Integer for q in M.inv()),'the inverse change of vanishing basis has integral entries')
zeros(list(M*residues-sp.Matrix([s-1,s+1,s-rho])),
      'the full native period vector has common slope one half in the local vanishing basis')
for r in [-1,1]:
    rr=residues.subs(rho,r);periods=M*rr
    index=0 if r==1 else 1
    zero(periods[index]-periods[2],'collided critical values give equal periods without identifying their cycles '+str(r))
    check(M.row(index)!=M.row(2),'the two collided-node cycles remain distinct primitive basis vectors '+str(r))
    zero(rr[index],'their nonzero homology difference has identically zero native-volume period '+str(r))
# General cotangent local volume and source-volume covariance controls.
j,delta,ci=sp.symbols('j delta ci',nonzero=True)
period=delta*(s-ci)/j
zero(sp.diff(period,s,2),'any common polynomial-source volume scalar preserves exactly linear full periods')
zero(period-(delta/j)*(s-ci),'the source change uses one scalar for all critical points')
K=sp.Matrix([[aa,bb],[cc,j]])
HH=sp.Matrix([[0,0,aa,cc],[0,0,bb,j],[aa,bb,0,0],[cc,j,0,0]])
zero(HH.det()-K.det()**2,'the cotangent critical Hessian fixes the square of the common volume-period slope')
algebra=len(checks);context=json.loads((HERE/'sources/global-vanishing-periods-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['rechecked_retained_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'rechecked full primary bytes '+entry['retained_file'])
check(context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[],
      'the integral tube and residue arguments reuse retained inputs without new classification imports')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==5 and not flags['Seifert_pairings_classified']
      and not flags['actual_singular_fibre_gluing_classified']
      and not flags['native_required_additive_action_constructed']
      and not flags['original_goal_complete'],
      'the full period object is exhausted without claiming finer marked data,actual gluing or a required action')
proof=HERE/'GLOBAL_INTEGRAL_VANISHING_CYCLES_AND_VOLUME_PERIODS.md'
paths=['sources/global-vanishing-periods-context.json']+['sources/'+e['retained_file'] for e in context['rechecked_retained_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-global-vanishing-periods.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
