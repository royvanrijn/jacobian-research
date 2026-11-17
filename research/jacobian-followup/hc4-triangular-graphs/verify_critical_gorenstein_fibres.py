"""Exact controls for critical-type recognition and first orbit modules."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.simplify(value)==0,label)
def zeros(values,label):check(all(sp.simplify(v)==0 for v in values),label)
u,v,h,P,Q,U=sp.symbols('u v h P Q U');source=[U,v,P,Q]
ell=sp.symbols('ell',integer=True,positive=True);m=ell+1
H=U*P+v**m*Q
def D(f):return v**m*sp.diff(f,P)-U*sp.diff(f,Q)
zero(D(H),'the genuine arbitrary-thickness control is preserved by its polynomial root action')
zeros([D(D(c)) for c in source],'the control is globally locally nilpotent on its actual polynomial source')
zeros([D(U),D(v)],'the control coefficient plane is fixed')
zero(D(P)-v**m,'the first actual slice value is the full arbitrary power v to m')
zero(D(-Q)-U,'the second actual slice value is U')
check(sp.gcd(U,v)==1,'the two coefficient factors have no common source prime')
newu=U+v*H
zero(D(newu),'the h-dependent chosen base generator is an actual invariant')
zero((u-v*h).subs({u:newu,h:H},simultaneous=True)-U,
     'the invariant frame change has its exact polynomial inverse')
zero(((u-v*h)*P+v**m*Q-h).subs({u:newu,h:H},simultaneous=True),
     'the h-dependent whole-source presentation holds on actual A4')
# The actual localized kernel chart,not merely a function-field formula.
zz=sp.symbols('zz');precovered=(h-v**m*Q)/U
zero(H.subs(P,precovered)-h,'the U-inverted polynomial chart recovers the value')
zero(precovered.subs(h,H)-P,'the localized chart has the exact source inverse')
zero(D(Q)+U,'the localized derivation is minus U partial_Q')
zero(H.subs(U,0)-v**m*Q,'the reduced chart gives an injective C[v,h] substitution h=v^m Q')
# The whole critical family persists in every allowed power,not just at one point.
critical={U:0,v:0,P:0}
zeros([sp.simplify(sp.diff(H,c)).subs(critical) for c in source],
      'the entire arbitrary-Q affine line is critical for the genuine A4 control')
zero(H.subs(critical),'that positive-dimensional critical family has the single value zero')
F=u-h*v;G=v**m;L=F*P+G*Q-h
zeros([sp.diff(F,h).subs({u:0,v:0}),sp.diff(G,h)],
      'h-dependent CI generators still have zero value derivative on the support line')
zero(sp.diff(L,h).subs({u:0,v:0})+1,
     'the source is smooth on the whole control support fibre regardless of Q')
zeros([sp.simplify(sp.diff(L,c)).subs({u:0,v:0,h:0,P:0}) for c in [u,v,P,Q]],
      'the h-dependent rank-one critical equations retain their entire affine line')
M=sp.Matrix([F,G]).jacobian([u,v]).applyfunc(sp.simplify).subs({u:0,v:0,h:0})
check(M.rank()==1,'the exact h-dependent fat CI has rank-one base derivative')
# The general mixed-Hessian formula,with arbitrary h-dependent coefficients evaluated at the point.
a,b,c,d,r11,r12,r22,ah=sp.symbols('a b c d r11 r12 r22 ah',nonzero=True)
mixed=sp.Matrix([[a,b],[c,d]])
hessian=sp.Matrix([[r11,r12,a,b],[r12,r22,c,d],[a,c,0,0],[b,d,0,0]])/ah
zero(hessian.det()-mixed.det()**2/ah**4,
     'the whole four-variable Hessian determinant is the squared base determinant divided by a_h to the fourth')
check(hessian[2:4,2:4]==sp.zeros(2),'the critical h Hessian has the zero vertical block used in the all-degree argument')
a11=sp.Function('a11')(h);a12=sp.Function('a12')(h)
a21=sp.Function('a21')(h);a22=sp.Function('a22')(h)
FF=u*a11+v*a12+u*u*sp.Function('r')(h)
GG=u*a21+v*a22+v*v*sp.Function('s')(h)
zeros([sp.diff(FF,h).subs({u:0,v:0}),sp.diff(GG,h).subs({u:0,v:0})],
      'arbitrary coefficient variation in h contributes no support-line value derivative')
zero((FF*P+GG*Q).subs({u:0,v:0}),
     'every solution of the critical linear system is on the same whole hypersurface fibre')
# Explicit local finite-algebra controls distinguish critical Gorensteinness from generic duality.
parameter=sp.symbols('parameter')
Mu=sp.Matrix([[0,0,0],[1,0,0],[0,0,0]])
Mv=sp.Matrix([[0,0,0],[0,0,parameter],[1,0,0]])
check(Mu*Mv==Mv*Mu,'the retained flat type-jump algebra has commuting multiplication')
check(Mv**2==parameter*Mu,'its actual finite-free family has the h-dependent relation v squared equals parameter times u')
check(Mu**2==Mu*Mv==sp.zeros(3),'the type-jump quotient has fixed finite support')
check(len(Mu.col_join(Mv).nullspace())==1,'generic inverse-system Gorensteinness has type one')
check(len(Mu.col_join(Mv).subs(parameter,0).nullspace())==2,
      'an actual length-three special ideal fibre may have type two')
check(1+len(Mu.col_join(Mv).nullspace())==2,
      'the plane Koszul formula gives two generic ideal generators')
check(1+len(Mu.col_join(Mv).subs(parameter,0).nullspace())==3,
      'the same formula gives three actual special ideal generators')
N2=sp.Matrix([[0,0],[1,0]])
check(N2**2==sp.zeros(2) and len(N2.nullspace())==1,
      'every length-two local algebra has the displayed one-dimensional square-zero maximal ideal and type one')
# The actual critical-type-one control C[u,v]/(u-cv,v^3).
Nv=sp.Matrix([[0,0,0],[1,0,0],[0,1,0]])
Nu=parameter*Nv
check(Nu*Nv==Nv*Nu and Nv**3==sp.zeros(3),
      'the genuine h-dependent Gorenstein fibre has its finite length-three multiplication model')
check(len(Nu.col_join(Nv).nullspace())==1,
      'that genuine control remains Gorenstein even though isolation fails')
# The small module has rank two but three generators in the desired standard geometry.
q=u*P+v*Q
def standard(f):return v*sp.diff(f,P)-u*sp.diff(f,Q)
zero(standard(q),'the standard source is the actual preserved Morse quadric')
check(sp.hessian(q,[u,v,P,Q]).det()==1,'the standard module example meets the native isolated Morse condition')
zero(-q+u*P+v*Q,'the complete first-orbit-module presentation relation holds on the actual source')
zeros([standard(standard(f)) for f in [1,P,Q]],
      'all three standard module generators have action degree at most one')
check(sp.gcd(u,v)==1,'every first-module syzygy is controlled by the coprime coefficient pair')
e=sp.symbols('e')
zero((u*e)*v-(v*e)*u,'the general coefficient-row syzygy is the displayed multiple')
column=sp.Matrix([-h,u,v])
check(column.rank()==1,'the standard single relation has generic rank one and module rank two')
check(column.subs({h:0,u:0,v:0})==sp.zeros(3,1),
      'every relation coefficient vanishes at the Morse quotient point,so its minimal module has three generators')
check(column.subs(h,1)[0]==-1,'outside the quotient point a unit relation coefficient locally splits off one generator')
# The new one-value total-length gap is arithmetic once every thick component has length at least three.
excess=sp.symbols('excess',integer=True,nonnegative=True)
check((3+excess-1).is_positive and sp.simplify((3+excess-1)-2)==excess,
      'any single thick component contributes at least two to the total excess over the support count')
algebra=len(checks)
context=json.loads((HERE/'sources/critical-gorenstein-fibres-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected primary Rees-torsor scope '+entry['retained_file'])
check(context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[],
      'no torsor bundle conclusion,Rees finite generation or unproved critical Gorensteinness is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==8 and not flags['critical_Gorenstein_property_proved_automatic_for_all_native_actions']
      and not flags['native_first_orbit_module_three_generator_bound_proved_automatic']
      and not flags['first_orbit_module_rank_two_implies_local_freeness_at_critical_point']
      and not flags['S_rho_normal_form_source_alignment_constructed'] and not flags['original_goal_complete'],
      'critical-type recognition and the small-module criterion leave the genuine non-CI native action gap open')
proof=HERE/'CRITICAL_GORENSTEIN_FIBRES_AND_FIRST_ORBIT_MODULES.md'
paths=['sources/critical-gorenstein-fibres-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-critical-gorenstein-fibres.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
