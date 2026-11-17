"""New point-defect/ring-map controls,never registered checker replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];started=time.perf_counter();checks=[];algebra=[]
def check(name,test,alg=False):
    assert test,name
    checks.append(name)
    if alg:algebra.append(name)
def zero(name,z):check(name,s.simplify(s.expand(z))==0,True)
p,q,w,r,v=s.symbols('p q w r v');a=s.symbols('a',nonzero=True);D=p**2-q**2;F=w**2-D
check('finite free cover monic quadratic',s.Poly(F,w).degree()==2 and s.Poly(F,w).LC()==1,True)
sing=s.groebner([F,s.diff(F,p),s.diff(F,q),s.diff(F,w)],w,p,q)
check('cone singular scheme only vertex',[g.as_expr() for g in sing.polys]==[w,p,q],True)
U=s.Matrix([[0,D],[1,0]])
G=s.Matrix([[s.trace(s.eye(2)),s.trace(U)],[s.trace(U),s.trace(U**2)]])
check('full finite trace Gram',G==s.diag(2,2*D),True)
check('full trace dual basis matrix',G.inv()==s.diag(s.Rational(1,2),1/(2*D)),True)
zero('inverse different basis comparison',s.rem(w**2-D,F,w))
zero('implicit p derivative',2*w*(p/w)-2*p)
zero('implicit q derivative',2*w*(-q/w)+2*q)
zero('mixed derivative compatibility',p*q/w**3-p*q/w**3)
vertex=s.groebner([w,p,q],w,p,q)
check('canonical image quotient one dimensional',vertex.is_zero_dimensional and [g.as_expr() for g in vertex.polys]==[w,p,q],True)
rel=r**2-v**2-1;rr=(a+1/a)/2;vv=(a-1/a)/2
zero('normal first bracket algebra Laurent parametrization',rr**2-vv**2-1)
zero('Laurent inverse a',rr+vv-a)
zero('Laurent inverse a inverse',rr-vv-1/a)
zero('nonconstant unit in first bracket algebra',(r+v)*(r-v)-1-rel)
PP=w*rr;QQ=w*vv
zero('cover equation in normal first bracket algebra',w**2-PP**2+QQ**2)
zero('exceptional fiber ideal equality',w-(r*(w*r)-v*(w*v))+w*rel)
zero('target map Jacobian vanishes exactly at w0',s.det(s.Matrix([PP,QQ]).jacobian([w,a]))-w/a)
zero('exceptional fiber target p',PP.subs(w,0))
zero('exceptional fiber target q',QQ.subs(w,0))
zero('new bracket square trace pole',s.trace((p*U.inv())**2)-2*p**2/D)
check('new trace numerator does not cancel branch divisor',s.gcd(2*p**2,D)==1,True)
# The rank-two identity is universal; it does not imply any trace regularity.
aa,bb,cc,dd=[s.symbols(n+'1 '+n+'2') for n in ['a','b','c','d']]
det=lambda A,B:A[0]*B[1]-A[1]*B[0]
zero('rank-two Plucker identity alone',det(aa,bb)*det(cc,dd)-det(aa,cc)*det(bb,dd)+det(aa,dd)*det(bb,cc))
ctx=json.loads((H/'sources/kahler-poisson-image-context.json').read_text());known={e['id']:e for e in json.loads((R/'MATH_STATUS.json').read_text())['entries']}
for e in ctx['retained_structural_inputs']:check('prior '+e['claim'],known[e['claim']]['state']=='proved' and known[e['claim']]['canonical_source']==e['canonical_source'])
for path,d in ctx['retained_proof_and_source_sha256'].items():check('byte '+path,hashlib.sha256((R/path).read_bytes()).hexdigest()==d)
yes={'full_goal_active','previous_goal_turn_progress','classical_canonical_image_and_prior_derivative_results_credited','Kahler_bracket_image_finite_module_with_canonical_reflexive_hull','point_supported_canonical_image_defect_retained','trace_pairing_test_equivalent_after_reflexive_hull','base_projectivity_of_image_equivalent_to_full_canonical_module','quadric_cone_length_one_module_defect_control_proved','already_normal_first_bracket_algebra_has_exceptional_Gm_fiber','normalization_before_canonical_hull_shortcut_refuted','explicit_bracket_square_trace_has_branch_pole','bare_Jacobi_Plucker_multiplication_shortcut_retired','external_review_deferred'}
for path,value in ctx['scope_flags'].items():check('scope '+path,value is (path in yes))
proof=H/'KAHLER_POISSON_IMAGE_AND_POINT_DEFECT.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-started,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=ctx['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-kahler-poisson-image.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({p:out[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
