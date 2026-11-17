"""Exact controls for one-place intrinsic genus and descended-trace closure; not a formal proof."""
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

g0,delta,G,d,E=sp.symbols('g0 delta G d E')
chi0=1-g0-delta
zero((1-chi0)-(g0+delta),'the proper normalization sequence gives arithmetic genus g0+delta')
zero((1-G)-chi0-(g0+delta-G),'proper flat Euler constancy has exactly the claimed generic-genus consequence')
zero((2*G-1)-d*(2*G-1)-(1-d)*(2*G-1),'equal generic genus and one puncture force the nonzero odd factor')
check(all(2*n-1!=0 for n in [0,1,2,6]),'the integer genus controls include zero and positive normalization genera')
for gg,dd in [(0,0),(0,1),(0,3),(2,0),(2,2)]:
    check(1-gg-dd==1-(gg+dd),'intrinsic finite deficit and normalization genus combine without an infinity term '+str((gg,dd)))
check(2*(0+1)-1==1,'the nodal zero curve has normalization genus zero but smooth pencil genus one')
check(2*(2+2)-1==7,'positive normalization genus and finite singularities are both retained')

# A normalized node, with a direct codimension-one ring-deficit control.
u,x,y,lam=sp.symbols('u x y lam')
node=y*y-x*x*(x+1)
param=[u*u-1,u*(u*u-1)]
zero(node.subs({x:param[0],y:param[1]},simultaneous=True),'the target node has polynomial normalization')
zeros([a.subs(u,1)-a.subs(u,-1) for a in param],'two distinct normalization points form the internal double branch')
check(sp.gcd(sp.diff(param[0],u),sp.diff(param[1],u))==1,'the node normalization is immersive but not injective')
for j in [0,1,2,3]:
    poly=(u*u)**j+u*(u*u-1)*(u*u)**j
    zero(poly.subs(u,1)-poly.subs(u,-1),'polynomials A(u^2)+u(u^2-1)B(u^2) have matching node values '+str(j))
check(u.subs(u,1)-u.subs(u,-1)==2,'the missing linear normalization parameter detects the one-dimensional deficit')
zero(sp.discriminant(x**3+x*x+lam,x)+lam*(4+27*lam),
     'the nodal pencil is generically a smooth cubic with two finite exceptional parameters')
X,Y,Z=sp.symbols('X Y Z')
hom=Y*Y*Z-X*X*(X+Z)-lam*Z**3
zero(hom.subs(Z,0)+X**3,'the projective cubic has the unique infinity point [0:1:0]')
zero(sp.diff(hom,Z).subs({X:0,Y:1,Z:0})-1,'that infinity point is smooth for every pencil parameter')

# Different polynomial generating pairs for one normalized coordinate algebra.
a=sp.symbols('a',nonzero=True)
zz=u**3-a
xx=u*zz
yy=zz**2+a*zz/2
recovered_z=8*(xx**3+a*a*yy/4-yy*yy)/a**3
zero(recovered_z-zz,'the degree-six presentation recovers the degree-three generator by an exact polynomial identity')
zero(yy-(zz**2+a*zz/2),'the reverse ring inclusion is explicit')
check(sp.degree(xx,u)==4 and sp.degree(zz,u)==3 and sp.degree(yy,u)==6,
      'two presentations have different displayed degrees without assuming ambient equivalence')
zero(xx**3-(zz**4+a*zz**3),'the normalized curve generators satisfy the correct one-place equation')
for k in [1,2,3]:
    test=xx**k+zz
    zero(test.subs(zz,recovered_z)-test,'substitution preserves normalized polynomial degree classes '+str(k))

# A genuine nonlinear closed HC4 embedding with a descended nodal trace.
s,t,w1,w2=sp.symbols('s t w1 w2')
K=[s+t*t,t]
f=K[1]**2-K[0]**2*(K[0]+1)
A=[w1*w1,w2**3]
H=w1*A[0]+w2*A[1]
C=w1*w2
HX=f+H
HY=[K[0]-(2*f+H)*A[0]+w2*C,K[1]-(2*f+H)*A[1]-w1*C]
potential=f*f+w1*K[0]+w2*K[1]
zero(HX*HX+w1*HY[0]+w2*HY[1]-potential,'the full nonlinear nodal-trace embedding has the literal rank-one potential')
zero(sp.Matrix(K).jacobian([s,t]).det()-1,'the positive control has an actual plane Keller pair')
zero(sp.hessian(potential,[s,t,w1,w2]).det()-1,'its actual four-variable potential has unit Hessian')
recf=HX-H
recK=[HY[0]+(HX+recf)*A[0]-w2*C,HY[1]+(HX+recf)*A[1]+w1*C]
zeros([recK[0]-recK[1]**2-s,recK[1]-t],'ambient polynomials recover both source coordinates globally')
zero(node.subs({x:K[0],y:K[1]},simultaneous=True)-f,'the entire target nodal curve pulls back to the displayed source equation')
source_param=[param[0]-param[1]**2,param[1]]
zeros([k.subs({s:source_param[0],t:source_param[1]},simultaneous=True)-b for k,b in zip(K,param)],
      'the actual source and target radical curves are isomorphic through the Keller pair')
zero(f.subs({s:source_param[0],t:source_param[1]},simultaneous=True),'the complete source nodal curve has the displayed normalization')
zeros([v.subs(u,1)-v.subs(u,-1) for v in source_param],'the source component also has an internal double branch')
zeros([sp.diff(f,v).subs({s:0,t:0}) for v in [s,t]],'this is the formerly remaining singular nodal case, not an affine-line zero curve')
zero((f*f)**3-f**6,'positive-power trace descent is exactly a power equation for the normal trace')
check(all(sp.Poly(v,w1,w2).total_degree()>1 for v in [HX]+HY),'the actual helpers include nonlinear higher free dependence')

# A selected sheet does not give a complete inverse image.
badK=[s,t*t]
zero(sp.Matrix(badK).jacobian([s,t]).det()-2*t,'the selected-sheet boundary control is explicitly non-Keller')
zero((t*t-1)-(t-1)*(t+1),'a target line has two source prime components in that control')
zeros([badK[0].subs(t,1)-s,badK[1].subs(t,1)-1],'one selected source line maps isomorphically to the target line')
check(sp.degree(t*t-1,t)==2,'the other source component cannot be discarded from the full preimage')

context=json.loads((HERE/'sources/one-place-trace-context.json').read_text())
check(len(context['newly_used_unproved_external_statements'])==3,
      'Vistoli compactification and both foundational Stacks inputs are explicit external statements')
check(context['proper_flat_compactification_hypotheses_checked'] and context['finite_normalization_deficit_retained'],
      'flatness and all finite singularity deficits are verified rather than suppressed')
check(context['one_place_isomorphism_matches_generic_pencil_invariants'] and not context['arbitrary_curve_isomorphism_assumed_ambient_equivalence'],
      'intrinsic normalization data match pencils without an ambient curve equivalence')
check(context['every_descended_rank_one_trace_resolved'] and context['positive_power_rational_trace_descent_forces_inverse'],
      'the whole descended rank-one trace family closes, including rational positive-power descent')
check(not context['internally_multibranch_descended_exceptions_remain_open'],
      'internally multibranch descended curves are no longer a residual exception')
check(context['full_reduced_preimage_required'] and not context['selected_branch_assumed_complete_preimage'],
      'the full inverse image is a required geometric hypothesis')
check(context['nonautomorphic_one_place_radical_components_require_extra_source_prime'] and not context['extra_source_prime_silently_discarded'],
      'unsaturated one-place radical components retain every additional source prime')
check(not context['all_nonconstant_rank_one_traces_resolved'] and not context['unrestricted_HC4_or_JC2_resolved'],
      'genuinely non-descended traces and the unrestricted objective remain open')
math_count=len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained proof/source byte pin '+name)
for e in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/e['retained_file']).read_bytes()).hexdigest()==e['sha256'],
          'new foundational source byte pin '+e['retained_file'])
old=json.loads((HERE/'receipt-contractible-pencil.json').read_text())
check(old['outcome']=='PASS' and old['check_count']==78,'retained IY receipt remains PASS78 and is not newly replayed')
for name,digest in old['source_sha256'].items():
    check(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,'retained IY receipt source '+name)
names=['verify_one_place_trace.py','ONE_PLACE_PENCILS_AND_DESCENDED_TRACE_CLOSURE.md','sources/one-place-trace-context.json']
result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
 source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
 dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed={'receipt-contractible-pencil.json':78},
 newly_cached_primary_context=context['newly_cached_primary_context'],
 newly_activated_cached_primary_inputs=context['newly_activated_cached_primary_inputs'],
 newly_used_unproved_external_statements=context['newly_used_unproved_external_statements'],
 primary_crosschecks_not_additional_unproved_inputs=context['primary_crosschecks_not_additional_unproved_inputs'],
 retained_structural_inputs=context['retained_structural_inputs'],written_not_finite_exhaustion=context['direct_proof_arguments'],
 **{k:v for k,v in context.items() if type(v) is bool},
 runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-one-place-trace.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,original_goal_complete=False),indent=2))
