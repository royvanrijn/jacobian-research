"""Exact algebra controls for the all-degree three-point module barrier."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
x,y,z,w=sp.symbols('x y z w');coords=[x,y,z,w];q=x*y+z*w;k=y*z-q*w
def deriv(values):return lambda f:sp.expand(sum(a*sp.diff(f,c) for a,c in zip(values,coords)))
D1=deriv([z,0,0,-y]);D2=deriv([w,0,-y,0]);D=deriv([y*z+q*w,0,-q*y,-y*y])
zeros([D1(D2(c))-D2(D1(c)) for c in coords],'the two actual root actions commute on all generators')
zeros([D(c)-y*D1(c)-q*D2(c) for c in coords],'the displayed action is the invariant combination yD1+qD2')
zeros([D1(y),D1(q),D2(y),D2(q)],'both coefficient functions belong to the two common kernels')
zeros([D(q),D(k)],'q and k are actual polynomial invariants')
zeros([D(D(x))+2*q*y*y,D(D(D(x))),D(D(z)),D(D(w))],
      'the explicit generator iterates terminate with their exact coefficients')
check(sp.hessian(q,coords).det()==1,'the actual quadric has unit Hessian determinant')
check(sp.Matrix([sp.diff(q,c) for c in coords]).jacobian(coords).rank()==4,
      'the quadric gradient has its sole zero at the origin')
zeros([sp.diff(D(c),a).subs(dict.fromkeys(coords,0)) for c in coords for a in coords],
      'the actual action has zero linearization at its Morse point')
check(sp.gcd(q,y)==1,'q and y have no common source prime')
zero(D(x).subs(y,0)-z*w*w,'the remaining action coefficient is not divisible by y')
check(sp.gcd_list([D(c) for c in coords])==1,'the actual quadric action is primitive')
Q,K=sp.symbols('Q K');zz=(K+Q*w)/y;xx=(Q-zz*w)/y
zero(q.subs({x:xx,z:zz},simultaneous=True)-Q,'the localized source chart recovers q')
zero(k.subs({x:xx,z:zz},simultaneous=True)-K,'the localized source chart recovers k')
zero(zz.subs({Q:q,K:k},simultaneous=True)-z,'the chart inverse recovers z')
zero(xx.subs({Q:q,K:k},simultaneous=True)-x,'the chart inverse recovers x')
zero(D(w)+y*y,'the localized action is minus y squared partial_w')
zero(sp.Matrix([q.subs(y,0),k.subs(y,0)]).jacobian([z,w]).det()+z*w*w,
     'modulo y the two remaining invariant generators are algebraically independent')
check(sp.Matrix([y,k]).jacobian(coords).subs(dict.fromkeys(coords,0)).rank()==1,
      'the proposed base projection has rank one at the Morse point')
zero(D(-w)-y*y,'the first exact plinth slice is y squared')
zero(D(-q*x-z*z)-q*k,'the second exact plinth slice is q times k')
zero(D(-z/q)-y,'after inverting the value the plinth contains y')
zero(q.subs(y,0)-z*w,'the quotient on y zero has q equal zw')
zero(k.subs(y,0)+q.subs(y,0)*w,'the quotient on y zero has k equal minus qw')
zero(k.subs({y:0,z:-K/w**2},simultaneous=True)-K,
     'every nonzero k fibre of the base has the displayed Laurent source chart')
r=sp.symbols('r')
intersection=sp.groebner([r*y,r*K,(1-r)*y,(1-r)*Q],r,y,Q,K)
check([p.as_expr() for p in intersection.polys if not p.as_expr().has(r)]==[y,Q*K],
      'the two distinct quotient lines intersect as the radical ideal (y,qk)')
# The actual three-point whole curve and both collision cases.
P,V,u,v,rho,t,s,alpha,beta=sp.symbols('P V u v rho t s alpha beta')
basecoords=[P,V,u,v];T=u*v+rho;B=u*(T*T-1)+T;h=2*P*V+B
Bu=sp.diff(B,u);Bv=sp.diff(B,v)
def dd(values,f):return sp.expand(sum(a*sp.diff(f,c) for a,c in zip(values,basecoords)))
du=[0,-Bu,2*P,0];dv=[0,-Bv,0,2*P]
zero(dd(du,h),'the first standard module generator fixes the actual potential')
zero(dd(dv,h),'the second standard module generator fixes the actual potential')
zeros([dd(du,dd(dv,c))-dd(dv,dd(du,c)) for c in basecoords],
      'the two actual suspension module generators commute')
zeros([dd(du,P),dd(dv,P)],'both standard actions fix the other coefficient variable')
d=t*t-1;uc=(s-t)/d;vc=(t-rho)*d/(s-t)
sub={u:uc,v:vc};E=2*s*t-t*t-1
A=sp.expand(d*(s-t)+(t-rho)*E)
zero(T.subs(sub,simultaneous=True)-t,'the rational whole-curve chart recovers its parameter')
zero(B.subs(sub,simultaneous=True)-s,'the rational curve is the exact level B=s')
zero(A-(-2*t**3+(3*s+rho)*t*t-2*s*rho*t+rho-s),'the exact cubic numerator is expanded correctly')
zero(Bu.subs(sub,simultaneous=True)-A/(s-t),'the u directional derivative has the full cubic numerator')
zero(Bv.subs(sub,simultaneous=True)-(s-t)*E/d**2,'the v directional derivative has the full quadratic numerator')
F=sp.expand(alpha*A*d*d+beta*(s-t)**2*E)
zero((alpha*Bu+beta*Bv).subs(sub,simultaneous=True)-F/((s-t)*d*d),
     'the fixed-line equation has exactly the stated nonvanishing denominator')
check(sp.degree(F,t)==7 and sp.Poly(F,t).LC()==-2*alpha,
      'the mixed direction has degree seven regardless of all coefficient degrees')
check(sp.degree((s-t)**2*E,t)==4,'the second direction has lower degree four for the primitive normalization argument')
zero(F.subs(t,s)-alpha*(s-rho)*(s*s-1)**3,'the mixed direction never vanishes at the moving puncture')
for r in [1,-1]:
    zero(F.subs(t,r)-beta*(s-r)**2*E.subs(t,r),'the mixed direction at the fixed puncture '+str(r))
zero(E.subs(t,1)-2*(s-1),'the pure v direction avoids the positive puncture')
zero(E.subs(t,-1)+2*(s+1),'the pure v direction avoids the negative puncture')
zero(E.subs(t,s)-(s*s-1),'the pure v direction avoids the moving puncture')
check(sp.degree(E,t)==2,'the pure v direction has a finite complex root')
zero(A.subs(t,s)-(s-rho)*(s*s-1),'the pure u direction avoids the moving puncture')
zero(A.subs(t,1)-2*(s-1)*(1-rho),'the pure u direction has the exact positive collision factor')
zero(A.subs(t,-1)-2*(s+1)*(1+rho),'the pure u direction has the exact negative collision factor')
check(sp.degree(A,t)==3 and sp.Poly(A,t).LC()==-2,'the pure u numerator always has degree three')
Qplus=-2*t*t+(3*s-1)*t+s-1;Qminus=-2*t*t+(3*s+1)*t-s-1
zero(A.subs(rho,1)-(t-1)*Qplus,'the rho=1 collision removes exactly the stated puncture factor')
zero(A.subs(rho,-1)-(t+1)*Qminus,'the rho=-1 collision removes exactly the stated puncture factor')
for label,poly,expected in [('positive',Qplus,[4*(s-1),-2*(s+1),s*s-1]),
                            ('negative',Qminus,[2*(s-1),-4*(s+1),s*s-1])]:
    check(sp.degree(poly,t)==2 and sp.Poly(poly,t).LC()==-2,'the '+label+' collision leaves a genuine quadratic')
    for r,value in zip([1,-1,s],expected):
        zero(poly.subs(t,r)-value,'the '+label+' collision residual avoids puncture '+str(r))
# Exact regular-fibre-free translation of the pure source alignment target.
X,Z,PP,QQ=sp.symbols('X Z PP QQ');targetcoords=[X,Z,PP,QQ]
C=(Z-1)*(Z-2)*(Z-3);a=(3*Z*Z-12*Z+10)/2
dr=((rho+3)*Z*Z-(3*rho+13)*Z+2*rho+12)/2;target=a*X*PP+C*QQ+dr
zero(sp.rem(a*sp.diff(C,Z)-1,C,Z),'the actual target coefficients have no common root')
zero(C*sp.diff(target,PP)-a*X*sp.diff(target,QQ),'the target coefficient translation fixes its potential')
for r,value in [(1,1),(2,-1),(3,rho)]:
    zero(dr.subs(Z,r)-value,'the target fixed plane has only critical value at root '+str(r))
    check(a.subs(Z,r)!=0,'the target action forces X zero at root '+str(r))
algebra=len(checks)
context=json.loads((HERE/'sources/three-point-translation-barrier-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected abstract-only primary context '+entry['retained_file'])
check(context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[],
      'no external classification or unproved hypothesis is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==6 and not flags['quadric_control_satisfies_fixed_cylindrical_plinth_support']
      and not flags['all_h_preserving_actions_or_actions_fixing_one_hyperbolic_variable_classified']
      and not flags['fixed_support_alone_proved_to_force_minimal_native_plinth_length']
      and not flags['S_rho_normal_form_source_alignment_constructed'] and not flags['original_goal_complete'],
      'entire standard quotient classes are excluded while the native classification,minimum length and alignment remain open')
proof=HERE/'THREE_POINT_TRANSLATION_QUOTIENT_BARRIER_AND_MORSE_COUNTERCONTROL.md'
paths=['sources/three-point-translation-barrier-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-three-point-translation-barrier.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
