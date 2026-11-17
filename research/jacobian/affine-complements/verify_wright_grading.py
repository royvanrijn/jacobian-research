"""Exact controls for all-height Wright grading and the finite trace boundary.

The universal proofs are written; plane torus linearization is imported
only for the invariant-target-algebra corollary.
"""
from pathlib import Path
from itertools import product
from time import monotonic
import sympy as S,json,hashlib,platform
HERE=Path(__file__).resolve().parent
start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(value)==0)
X,Y,Z,n,l,e=S.symbols('X Y Z n l e')
f,g,df,dg=S.symbols('f g df dg')
formal=(n*f+2*Z*df)*(l*g+Z*dg)-Z*df*(e*g+2*Z*dg)
zero('universal differentiated homogeneous-pair coefficient',
     formal-(n*l*f*g+n*Z*f*dg+(2*l-e)*Z*df*g))
T=X*Y;B=X*X*Y;D=X**3*Y+X
zero('graded relation T squared',T*T-Y*B)
zero('graded relation T times D',T*D-B*(B+1))
zero('graded relation Y times D',Y*D-T*(B+1))
zero('graded relation Y times D squared',Y*D*D-B*(B+1)**2)

# Group bounded monomial checks by weight; the all-height normal-form
# argument is separately written and is not inferred from this list.
controls={}
for i,j,k,d in product(range(6),repeat=4):
    if i+j+k+d>5:continue
    weight=2*i+j-d
    monomial=Y**i*T**j*B**k*D**d
    if weight<0:generator=D**(-weight)
    elif weight%2==0:generator=Y**(weight//2)
    else:generator=T*Y**((weight-1)//2)
    quotient=S.cancel((monomial/generator).subs(Y,Z/X**2))
    controls.setdefault(weight,[]).append(quotient)
for weight,values in sorted(controls.items()):
    check('full-piece monomial controls at weight '+str(weight),
          all(not value.has(X,Y) and S.denom(value)==1 and value.is_polynomial(Z) for value in values))
# Independent actual differentiation controls for both parities.
for integer in range(9):
    ell=integer//2+1;epsilon=1 if integer%2==0 else 0
    ff=(Z+1)**integer*(2+3*Z+5*Z**2);gg=7+11*Z
    PP=X**integer*ff.subs(Z,B);QQ=X**epsilon*Y**ell*gg.subs(Z,B)
    actual=S.diff(PP,X)*S.diff(QQ,Y)-S.diff(PP,Y)*S.diff(QQ,X)
    expected=Z**(ell-1)*(integer*ell*ff*gg+integer*Z*ff*S.diff(gg,Z)+(integer+1)*Z*S.diff(ff,Z)*gg)
    zero('actual homogeneous Jacobian parity control '+str(integer),actual-expected.subs(Z,B))
    check('weight sum in parity control '+str(integer),2*ell-epsilon==integer+1)
zero('the only ordinary homogeneous resonance uses forbidden X',
     S.diff(X,X)*S.diff(Y,Y)-S.diff(X,Y)*S.diff(Y,X)-1)
zero('regular negative weight generator loses constant Jacobian',
     S.diff(D,X)-(1+3*B))

# Pole-filtration raising bound via universal V1 identities.
delta=lambda H:S.expand(-X*S.diff(H,X)+2*Y*S.diff(H,Y))
C=1+B+D*D;AA=T
old=[-Y-C*C-2*AA*D-B*D*D,-AA-D-2*B*D-D**3,C,D]
new=[2*old[0]+2*C*C-4*old[1]*D,old[1]-2*D+4*C*D,-2*D*D,-D]
for j,(value,expected) in enumerate(zip(old,new)):zero('V1 grading derivative '+str(j),delta(value)-expected)
for value,weight in [(Y,2),(T,1),(B,0),(D,-1)]:zero('actual ring eigenweight '+str(weight),delta(value)-weight*value)
zero('actual source-volume divergence',S.diff(-X,X)+S.diff(2*Y,Y)-1)

saved=json.loads((HERE/'wright_grading_certificates.json').read_text())
check('producer source seal',saved['source_sha256']==hashlib.sha256(
    (HERE/'compute_wright_grading.py').read_bytes()).hexdigest())
counter=saved['trace_countercontrol']
x,w,p,q,t=S.symbols('x w p q t')
local={str(v):v for v in (X,Y,Z,x,w,p,q,t)}
parse=lambda value:S.sympify(value,locals=local)
P=parse(counter['P']);Q=parse(counter['Q'])
J=S.expand(S.diff(P,x)*S.diff(Q,w)-S.diff(P,w)*S.diff(Q,x))
zero('countercontrol is ramified along the single smooth line',J-w)
E=lambda H:S.expand((-x*S.diff(H,x)+2*w*S.diff(H,w))/3)
zero('countercontrol grading preserves ramification ideal',E(w)-S.Rational(2,3)*w)
zero('countercontrol volume is homogeneous of weight one',
     E(J)+(S.diff(-x/3,x)+S.diff(2*w/3,w))*J-J)
F=parse(counter['monic_quartic'])
zero('finite quartic relation on actual source',F.subs({p:P,q:Q},simultaneous=True))
zero('polynomial inverse recovering w',P.subs(w,p-x-x**3)-p)
zero('polynomial inverse recovering q',Q.subs(w,p-x-x**3)-q-F/4)
check('finite module has exact monic rank four',S.Poly(F,x).degree()==4 and S.Poly(F,x).LC()==1)
# Power traces from Newton sums, independent of the producer's matrix powers.
coeff=S.Poly(F,x)
newton={0:S.Integer(4)}
for j in range(1,5):
    newton[j]=-sum(coeff.nth(4-i)*newton[j-i] for i in range(1,j))-j*coeff.nth(4-j)
for j in range(5):zero('independent Newton trace '+str(j),newton[j]-parse(counter['power_traces'][str(j)]))
M=S.Matrix([[parse(value) for value in row] for row in counter['multiplication_matrix_x']])
zero('multiplication matrix characteristic polynomial',M.charpoly(x).as_expr()-F)
def trace_poly(H):
    pp=S.Poly(H,x,domain=S.QQ.frac_field(p,q))
    return S.cancel(sum(coef*newton[powers[0]] for powers,coef in pp.terms()))
Fp=S.Poly(F,x,domain=S.QQ.frac_field(p,q))
inv_der=S.invert(S.Poly(S.diff(F,x),x,domain=S.QQ.frac_field(p,q)),Fp).as_expr()
for parameter,numerator in [(p,4*x),(q,S.Integer(4))]:
    for j in range(1,5):
        lifted=S.rem(j*x**(j-1)*numerator*inv_der,F,x)
        zero('trace commutes with lifted partial '+str(parameter)+' '+str(j),
             trace_poly(lifted)-S.diff(newton[j],parameter))
DP=S.rem(E(P).subs(w,p-x-x**3),F,x)
DQ=S.rem(E(Q).subs(w,p-x-x**3),F,x)
TP=parse(counter['target_trace_P']);TQ=parse(counter['target_trace_Q'])
zero('independent field trace of grading P',trace_poly(DP)-TP)
zero('independent field trace of grading Q',trace_poly(DQ)-TQ)
zero('target trace divergence exactly four',S.diff(TP,p)+S.diff(TQ,q)-4)
pc,qc=map(parse,counter['branch_parameterization'])
zero('branch parametrization P',P.subs({x:t,w:0})-pc)
zero('branch parametrization Q',Q.subs({x:t,w:0})-qc)
zero('branch normalization has residue degree one',parse(counter['branch_rational_inverse']).subs({p:pc,q:qc})-t)
transverse=TP.subs({p:pc,q:qc})*S.diff(qc,t)-TQ.subs({p:pc,q:qc})*S.diff(pc,t)
zero('complete trace transversality polynomial',transverse-parse(counter['branch_transversality']))
check('trace is transverse at the smooth branch origin',S.cancel(transverse.subs(t,0))==-8)
disc=parse(counter['discriminant'])
zero('exact quartic discriminant',S.discriminant(F,x)-disc)
check('branch discriminant vanishes at control origin',disc.subs({p:0,q:0})==0)
check('trace does not preserve the branch ideal',
      (TP*S.diff(disc,p)+TQ*S.diff(disc,q)).subs({p:0,q:0})!=0)
# At target zero, the two affine sheets x^2=-2 have w=x and each normal
# contribution equals4; do not discard them from the finite trace.
zero('two unramified-sheet contributions retained',
     S.rem(E(Q).subs(w,x)-4,x*x+2,x))
zero('branch-origin normal trace comes from affine sheets',TQ.subs({p:0,q:0})-2*4)

manifest_path=HERE/'sources/grading-linearization-manifest.json'
manifest=json.loads(manifest_path.read_text());source=manifest['source']
blob=(manifest_path.parent/source['file']).read_bytes()
check('linearization primary source byte pin',len(blob)==source['bytes'])
check('linearization primary source sha256',hashlib.sha256(blob).hexdigest()==source['sha256'])
check('core mate obstruction is independent of linearization',
      manifest['homogeneous_mate_obstruction_uses_linearization'] is False)
check('linearization proof remains external',manifest['external_theorem_independently_reproved'] is False)
check('correct primary Corollary4.6 retained',b'Corollary 4.6' in blob and b'linearizable' in blob)
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    bounded_weight_module_controls=sum(map(len,controls.values())),
    written_results=['no homogeneous element of the Wright ring has a Jacobian mate in the ring',
      'strict weight-range and commuting-extreme constraints for every Keller target presentation',
      'no grading-invariant polynomial Keller target algebra, using external plane torus linearization',
      'grading raises B-pole order by at most one',
      'finite-cover grading trace is polynomial with divergence four; branch tangency is not automatic'],
    trace_countercontrol=dict(finite_degree=4,ramification_prime='w=0',residue_degree=1,
       normalized_target_divergence=4,branch_tangency=False,transversality_at_origin='-8',
       is_wright_candidate=False,is_keller_counterexample=False),
    subsequent_all_height_pole_three_exclusion="PF2RULED8 (separate intrinsic proof and receipt)",
            frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
    assurance=dict(all_height_homogeneous_obstruction='written proof',
       invariant_target_algebra_obstruction='written proof plus external Gutwirth linearization',
       trace_divergence='written general field proof with exact quartic control',
       actual_wright_trace_tangency_proved=False,unrestricted_mixed_weight_pairs_excluded=False,
       finite_map_on_wright_model_constructed=False,actual_ruled_compatibility_proved=False,
       quartic_packets_excluded=False,formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
       Path(__file__),HERE/'compute_wright_grading.py',HERE/'wright_grading_certificates.json',
       HERE/'WRIGHT_GRADING_OBSTRUCTION.md',HERE/'TRACE_GRADING_BOUNDARY.md',manifest_path]},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'wright_grading_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; homogeneous mates excluded at all heights; trace branch tangency remains an extra gate')
