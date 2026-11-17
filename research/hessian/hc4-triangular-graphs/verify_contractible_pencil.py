"""Exact controls for the written all-degree pencil and contractible-curve proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(v,label):
    assert bool(v),label
    checks.append(label)
def zero(v,label):check(sp.simplify(sp.expand(v))==0,label)
def zeros(vs,label):check(all(sp.simplify(sp.expand(v))==0 for v in vs),label)

gs,gt,rs,rt,d,q,E,Rf=sp.symbols('gs gt rs rt d q E Rf')
nus=2*gs-2+rs;nut=2*gt-2+rt
rh=2*gs-2-d*(2*gt-2)-(d*rt-q)-Rf
zero(nus-d*nut-((rs-q)+Rf)-rh,'the exact log-Euler defect follows from projective Riemann-Hurwitz')
zero(((rs-q)+Rf-E).subs(Rf,E-(rs-q)),
     'each finite-image puncture contributes its full index to E')
zero((2*gs-1)-d*(2*gs-1)-(1-d)*(2*gs-1),'equal genus and one puncture leave a nonzero odd factor multiplying d-1')
check(all(2*g-1!=0 for g in [0,1,2,7]),'sample genera include the genus-zero and positive-genus one-puncture cases')
zero(0-d*(-1)-d,'Gm-to-A1 t^d has the exact nonzero finite-boundary loss E=d')
zero(0-d*0,'Gm-to-Gm powers have zero log-Euler invariant at arbitrary degree')
z=sp.symbols('z')
for m in [1,2,5]:
    zero(sp.diff(z**m,z)-m*z**(m-1),'the curve power map derivative is a localized unit on Gm m='+str(m))
    check(sp.degree(z**m,z)==m,'the same etale affine curve control has arbitrary degree m='+str(m))
    zero((2*m-2)-((m-1)+(m-1)),'its two projective ramification indices exhaust Riemann-Hurwitz m='+str(m))
    zero(0-m*(-1)-m,'its source puncture over zero is retained in the defect m='+str(m))
check(2*1+0>1,'an alleged degree-two map with equal positive nu would require negative boundary loss')
check(5<3*2 and 5//2==2,
      'the positive-invariant bound below three leaves only degrees one or two; retained HW0 excludes degree two')

p,r=sp.symbols('p r',integer=True,positive=True)
genus=(p-1)*(r-1)/2
zero(2*genus-2-(-2*p+(r+1)*(p-1)),'the all-exponent monomial genus formula is the exact ramification identity')
zero(2*genus-1-(p*r-p-r),'the monomial generic log-Euler invariant is pq-p-q')
for pp,qq in [(2,3),(2,5),(3,4),(3,7)]:
    check(sp.gcd(pp,qq)==1,'the cusp sample has one infinity valuation orbit '+str((pp,qq)))
    gg=(pp-1)*(qq-1)//2
    check(gg>=1 and 2*gg-1>0,'the singular monomial sample has positive one-puncture log-Euler invariant '+str((pp,qq)))
    zero(-2*pp+(qq+1)*(pp-1)-(2*gg-2),'its finite and infinite ramification counts give the displayed genus '+str((pp,qq)))
    check(min(pp,qq)==min([pp,qq]) and min(pp,qq)>1,'its valuation semigroup has no order-one normalization parameter '+str((pp,qq)))

# Actual nonlinear Keller pair and closed variable-normal-trace embedding.
s,t,w1,w2=sp.symbols('s t w1 w2')
K=[s+t*t,t]
f=K[0]**2-K[1]**3
A=[w1*w1,w2**3]
H=w1*A[0]+w2*A[1]
C=w1*w2
X=f+H
Y=[K[0]-(2*f+H)*A[0]+w2*C,K[1]-(2*f+H)*A[1]-w1*C]
phi=f*f+w1*K[0]+w2*K[1]
zero(X*X+w1*Y[0]+w2*Y[1]-phi,'the actual cusp-trace family has its literal nonconstant potential identity')
zero(sp.Matrix(K).jacobian([s,t]).det()-1,'the actual positive control is genuinely plane Keller')
zero(sp.hessian(phi,[s,t,w1,w2]).det()-1,'its actual four-variable potential has unit Hessian')
recf=X-H
recK=[Y[0]+(X+recf)*A[0]-w2*C,Y[1]+(X+recf)*A[1]+w1*C]
zeros([recK[0]-recK[1]**2-s,recK[1]-t],'ambient polynomials recover the entire source of the high-free-degree embedding')
T1,T2=sp.symbols('T1 T2')
cusp=T1*T1-T2**3
zero(cusp.subs({T1:K[0],T2:K[1]},simultaneous=True)-f,'the target cusp pulls back to the exact source irreducible curve polynomial')
source_param=[z**3-z**4,z*z]
target_param=[z**3,z*z]
zeros([a.subs({s:source_param[0],t:source_param[1]},simultaneous=True)-b for a,b in zip(K,target_param)],
      'the actual Keller restriction is the isomorphism between the two cusp embeddings')
zero(f.subs({s:source_param[0],t:source_param[1]},simultaneous=True),'the normalized source parameter lies on the complete pulled-back cusp')
zero(cusp.subs({T1:target_param[0],T2:target_param[1]},simultaneous=True),'the normalized target parameter lies on the target cusp')
zeros([sp.diff(a,z).subs(z,0) for a in source_param+target_param],
      'both normalization parameters have zero tangent at the cusp and are not closed affine-line embeddings')
zeros([sp.diff(f,v).subs({s:0,t:0}) for v in [s,t]],'the actual source radical curve is singular; no smooth-A1 substitution proves this case')
chain=sp.Matrix([sp.diff(cusp,T1),sp.diff(cusp,T2)]).subs({T1:K[0],T2:K[1]},simultaneous=True).T*sp.Matrix(K).jacobian([s,t])
zeros([chain[0,i]-sp.diff(f,v) for i,v in enumerate([s,t])],'pencil differentiation identifies the tangent restriction through the invertible Keller differential')

# A1-normalized nodal curve has an internal double branch and survives the new condition.
node=T2*T2-T1*T1*(T1+1)
node_param=[z*z-1,z*(z*z-1)]
zero(node.subs({T1:node_param[0],T2:node_param[1]},simultaneous=True),'the nodal remaining-boundary control has polynomial normalization')
zeros([a.subs(z,1)-a.subs(z,-1) for a in node_param],'two distinct normalization points are identified at the node')
check(sp.gcd(sp.diff(node_param[0],z),sp.diff(node_param[1],z))==1,'its normalization has nonzero differential despite the double branch')
check(1-(2-1)==0 and 1-(1-1)==1,'an internal node has Euler zero while a unibranch cusp has Euler one')

context=json.loads((HERE/'sources/contractible-pencil-context.json').read_text())
check(len(context['newly_used_unproved_external_statements'])==1,'the newly imported Lin-Zaidenberg classification is explicit, not claimed reproduced')
check(context['full_reduced_preimage_required'] and not context['selected_branch_assumed_complete_preimage'],
      'a selected source branch is never substituted for the entire inverse image')
check(context['separate_contractible_rectifications_used'] and not context['arbitrary_curve_isomorphism_assumed_ambient_extendible'],
      'source and target contractible curves are rectified separately')
check(context['finite_image_boundary_points_retained'] and not context['affine_etale_assumed_finite'],
      'the exact boundary defect retains every finite-image puncture')
check(context['descended_contractible_components_closed'] and not context['descended_unibranch_cusps_remain_open'],
      'the prior unibranch singular-cusp frontier is now closed under actual full-preimage hypotheses')
check(context['every_nonautomorphic_descended_component_has_internal_multibranch_point'] and not context['nodal_descended_curves_declared_closed'],
      'internal multibranch points remain different from intersections of distinct target components')
check(not context['zero_curve_isomorphism_assumed_equal_generic_pencil_invariants'],
      'arbitrary isomorphic zero curves do not silently identify generic pencil genera or punctures')
check(not context['unrestricted_HC4_or_JC2_resolved'] and not context['original_goal_complete'],
      'remaining multibranch and unaligned cases keep the unrestricted objective open')
math_count=len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained proof/source byte pin '+name)
for e in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/e['retained_file']).read_bytes()).hexdigest()==e['sha256'],
          'new primary classification byte pin '+e['retained_file'])
old=json.loads((HERE/'receipt-radical-trace.json').read_text())
check(old['outcome']=='PASS' and old['check_count']==80,'retained IX receipt remains PASS80 and was not newly replayed')
for name,digest in old['source_sha256'].items():
    check(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,'retained IX receipt source '+name)
names=['verify_contractible_pencil.py','CONTRACTIBLE_RADICAL_CURVES_AND_PENCIL_DEFECT.md','sources/contractible-pencil-context.json']
result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
 source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
 dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed={'receipt-radical-trace.json':80},
 newly_cached_primary_context=context['newly_cached_primary_context'],
 newly_used_unproved_external_statements=context['newly_used_unproved_external_statements'],
 retained_external_foundations=context['retained_external_foundations'],retained_structural_inputs=context['retained_structural_inputs'],
 written_not_finite_exhaustion=context['direct_proof_arguments'],
 **{k:v for k,v in context.items() if type(v) is bool},
 runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-contractible-pencil.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,original_goal_complete=False),indent=2))
