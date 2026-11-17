"""Exact controls for the global torus gate and coordinate visibility proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
x,y,z,w,u,v,P,V,t,rho,alpha,a,b=sp.symbols('x y z w u v P V t rho alpha a b')
K2=w*(2*x*y+z*z-1)+z;vs=[x,y,z,w];H2=sp.factor(sp.hessian(K2,vs).det())
zero(H2-16*w*w*(2*x*y+z*z),'K2 full native Hessian determinant')
rw=[sp.diff(K2,q) for q in [x,y,z]]
zeros([rw[0]-2*w*y,rw[1]-2*w*x,rw[2]-(2*w*z+1)],'K2 w-relative ideal has the three stated exact generators')
zero(x*rw[2]-z*rw[1]-x,'the w-relative ideal contains x without localization')
zero(y*rw[2]-z*rw[0]-y,'the w-relative ideal contains y without localization')
zero((-2*z)*w-(1-(2*w*z+1)),'w is an actual unit on the whole relative torus')
torus={x:0,y:0,z:-1/(2*w)}
zeros([q.subs(torus,simultaneous=True) for q in rw],'the whole w-torus parameterization solves every relative equation')
zero(H2.subs(torus,simultaneous=True)-4,'the relative Hessian determinant passes identically on the K2 torus')
zero(K2.subs(torus,simultaneous=True)+w+1/(4*w),'the optional critical-value envelope agrees on the whole torus')
rz=[sp.diff(K2,q) for q in [x,y,w]]
zeros([rz[0]-2*w*y,rz[1]-2*w*x,rz[2]-(2*x*y+z*z-1)],'K2 z-relative scheme has its stated ideal')
zeros([rz[0].subs(w,0),rz[1].subs(w,0),rz[2].subs(w,0)-(2*x*y+z*z-1)],
      'the entire dimension-two quadric is contained in the z-relative scheme')
check(sp.Poly(2*x*y+z*z-1,x,y,z).total_degree()==2,'the nonzero quadric equation cuts a dimension-two locus in three variables')
for sign in [-1,1]:
    point={x:0,y:0,z:sign,w:-sp.Rational(sign,2)}
    zeros([sp.diff(K2,q).subs(point) for q in vs],'K2 full critical point '+str(sign))
    zero(H2.subs(point)-4,'K2 full Morse determinant '+str(sign))
h0=2*x*y+(z*z-1)*w;h1=h0+z
zeros([sp.diff(h0,q)-s for q,s in zip([x,y,z],[2*y,2*x,2*w*z])],
      'h0 w-relative scheme is the full singular union V(x,y,zw)')
check(sp.Matrix([2*y,2*x,2*w*z]).jacobian(vs).subs({x:0,y:0,z:0,w:0}).rank()==2,
      'h0 relative crossing has tangent dimension two instead of one')
zeros([sp.diff(h1,q)-s for q,s in zip([x,y,z],[2*y,2*x,2*w*z+1])],
      'h1 whole w-relative scheme is exactly the same torus')
for sign in [-1,1]:
    point={x:0,y:0,z:sign,w:-sp.Rational(sign,2)}
    zeros([sp.diff(h1,q).subs(point) for q in vs],'h1 full critical point '+str(sign))
    zero(sp.hessian(h1,vs).det().subs(point)-16,'h1 full Morse determinant '+str(sign))
# The two possible torus monomial maps have two points over every nonexceptional line value.
for exponent in [2,-2]:
    F=b+a*t**exponent
    zero(sp.diff(F,t)-a*exponent*t**(exponent-1),'Laurent monomial restriction has an everywhere nonzero differential on Gm,exponent='+str(exponent))
    numerator=sp.together(F-alpha).as_numer_denom()[0]
    check(sp.Poly(numerator,t).degree()==2,'the whole line restriction has generic degree two,exponent='+str(exponent))
    check(sp.discriminant(numerator,t)!=0,'the generic two points are reduced,exponent='+str(exponent))
T=u*v+rho;B=u*T*T-u+u*v+rho;S=2*P*V+B
Bu=sp.diff(B,u);Bv=sp.diff(B,v);HS=sp.factor(sp.hessian(S,[u,v,P,V]).det())
zero(Bu-(3*u*u*v*v+4*rho*u*v+v+rho*rho-1),'S_rho v-relative equation is the exact univariate-parameter family')
zero(Bv-u*(2*u*T+1),'S_rho u-relative scheme splits by comaximal factors')
zero((2*u*T+1)-2*T*u-1,'the line and torus factors have an explicit polynomial Bezout identity')
zero(Bv.subs(u,0),'the entire actual affine v-line is relatively critical retaining u')
vtorus=-1/(2*u*u)-rho/u
zero((2*u*T+1).subs(v,vtorus),'the second whole component has an actual Gm_u parameterization')
zero(u*(-2*T)-(1-(2*u*T+1)),'u is a regular unit on the second component')
Y=1-3*t*t+2*rho*t;curve={u:(t-rho)/Y,v:Y}
zero(Bu.subs(curve,simultaneous=True),'the entire parameter-uniform open curve lies on the v-relative scheme')
zero(T.subs(curve,simultaneous=True)-t,'the relative curve parameterization has its regular inverse t=uv+rho')
check(sp.Poly(Y,t).LC()==-3,'the retained curve open set is dense for every rho with no exceptional parameter')
H=-3*t**4+4*rho*t**3+6*t*t-12*rho*t+4*rho*rho+1
zero(HS.subs(curve,simultaneous=True)-4*H/Y**2,'the native Hessian class on the entire open curve is exactly 4H/Y^2')
check(sp.Poly(H,t).LC()==-3,'the relative determinant numerator has parameter-independent leading coefficient minus three')
check(sp.Poly(Y**2,t).LC()==9,'the would-be constant-four determinant denominator has leading coefficient nine')
check(sp.Poly(H-Y**2,t).LC()==-12,'the constant-four identity is impossible for every complex rho')
critical=[{u:0,v:1-rho*rho},{u:-sp.Rational(1,2),v:2*(rho-1)},
          {u:sp.Rational(1,2),v:-2*(rho+1)}]
for index,point in enumerate(critical):
    point.update({P:0,V:0})
    zeros([sp.diff(S,q).subs(point,simultaneous=True) for q in [u,v,P,V]],'S full critical point '+str(index))
    zero(HS.subs(point,simultaneous=True)-4,'S full Morse determinant is four for all rho at point '+str(index))
check(len({point[u] for point in critical})==3,'the three full critical points remain distinct even when critical values collide')
# Affine pullbacks are normalized only by a constant affine source map.
L=x+2*y+3*z-w+7;normal={x:x-2*y-3*z+w-7,y:y,z:z,w:w}
zero(L.subs(normal,simultaneous=True)-x,'an affine coordinate pullback has an actual constant affine normalizer')
zero(sp.Matrix([normal[q] for q in vs]).jacobian(vs).det()-1,'affine normalization preserves unit Hessian up to its constant determinant square')
algebra=len(checks);context=json.loads((HERE/'sources/torus-coordinate-visibility-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected primary bytes '+entry['retained_file'])
check(len(context['newly_imported_proved_external_statements'])==1 and context['newly_imported_proved_external_statements'][0]['hypotheses_checked'],
      'the sole new import has checked plane nonsingular hypotheses')
check(context['newly_used_unproved_external_statements']==[],'no unknown two-point fibre exclusion or full-source alignment is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==5 and not flags['arbitrary_two_point_Keller_fibre_excluded']
      and not flags['arbitrary_torus_component_substituted_for_whole_fibre']
      and not flags['native_K3_x_z_affine_pullbacks_excluded']
      and not flags['arbitrary_full_source_affine_alignment_proved'] and not flags['original_goal_complete'],
      'complete one-coordinate visibility excludes no unproved arbitrary-source or two-point Keller class')
proof=HERE/'TORUS_CRITICAL_CURVES_AND_ALL_COORDINATE_VISIBILITY.md'
source_paths=['sources/torus-coordinate-visibility-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in source_paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-torus-coordinate-visibility.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
