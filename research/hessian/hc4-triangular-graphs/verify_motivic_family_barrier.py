"""Exact chart and hypothesis controls;this does not implement motivic homotopy."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
x,y,z,w,u,v,s,t,rho,X,P,Q,F,G,d,C,aa,bb,ee=sp.symbols('x y z w u v s t rho X P Q F G d C aa bb ee')
# Actual entire cotangent family charts,for any zero-fibre cardinality and any d.
eq=P*F+Q*G+d-s
zero(eq.subs(P,(s-d-Q*G)/F),'the F-open chart solves the whole cotangent value equation')
zero(eq.subs(Q,(s-d-P*F)/G),'the G-open chart solves the whole cotangent value equation')
zero(P/G+Q/F-(eq+s-d)/(F*G),'the exact Cech slice difference retains the entire value-dependent class')
zero(G*sp.diff(eq,P)-F*sp.diff(eq,Q),'the global translation preserves the family equation')
# The universal N-point complement model uses a formal coefficient C_N(z).
model=X*P+C*Q-1
zero(model.subs(P,(1-C*Q)/X),'W_N X-open actual affine-line chart')
zero(model.subs(Q,(1-X*P)/C),'W_N C_N-open actual affine-line chart')
zero(model.subs({X:1,P:1,Q:0},simultaneous=True),'the global z-section works identically for every C_N and N')
zero(model.subs(C,0)-(X*P-1),'each simple root has the whole fibre Gm times A1')
zero((X*P-1).subs(P,1/X),'the exceptional Gm fibre has its exact inverse coordinate')
zero(model.subs(C,0)-P*sp.diff(model,P).subs(C,0)+1,'X and P cannot both vanish on an exceptional fibre')
check(sp.diff(model,Q)==C,'outside the roots the Q derivative gives the trivial smooth A2 bundle')
zero((X*P+C*Q).subs({P:s*P,Q:s*Q},simultaneous=True)-s*(X*P+C*Q),
     'rescaling one isotropic plane trivializes the nonzero quadratic value family')
zeros([(s*P)/s-P,(s*Q)/s-Q],'both inverse quadratic rescalings are regular when s is a unit')
# General D division by z^2-1;the remainder parameters are unrestricted.
Cz=z*z-1;D=aa*z+bb+Cz*ee;p=s-D
Delta=(s-bb)**2-aa**2
zero(Delta-(p+Cz*ee)*(s-bb+aa*z)-aa**2*Cz,
     'the universal Bezout identity makes z^2-1 a unit on the entire h_D exceptional curve over S')
zero(Delta-(s-(aa+bb))*(s-(-aa+bb)),'the only inverse denominators are the two excluded critical values')
B=Cz*w+D
zero((B-s).subs(w,(s-D)/Cz),'the whole h_D exceptional curve has its exact punctured-line formula')
zero((B-s).subs({z:0,w:bb-ee-s},simultaneous=True),'h_D has the regular t0=0 section for arbitrary D')
check(sp.diff(B,w)==Cz,'the h_D plane curve is relatively smooth because its w derivative is a unit')
K2=w*(2*x*y+z*z-1)+z;native_u=w;native_v=1-2*x*y-z*z
zero(K2+native_u*native_v-z,'the entire K2 family has exact reverse z=s+uv')
zz=s+u*v;Fs=v+zz*zz-1
zero(Fs.subs({s:K2,u:native_u,v:native_v},simultaneous=True)+2*x*y,
     'the entire K2 family satisfies exactly its uniform suspension relation')
zero(1-2*x*((1-zz*zz-v)/(2*x))-zz*zz-v,'the localized reverse y expression recovers v')
zero(Fs-v*(1+u*(zz+s))-s*s+1,'v is a genuine unit on the entire exceptional curve when s^2-1 is a unit')
kcurve={u:(z-s)/(1-z*z),v:1-z*z}
zero(Fs.subs(kcurve,simultaneous=True),'K2 whole exceptional plane curve parameterization')
zero(zz.subs(kcurve,simultaneous=True)-z,'K2 curve regular inverse is z=s+uv')
zeros([kcurve[u].subs(z,zz)-u-u*Fs/(1-zz*zz),kcurve[v].subs(z,zz)-v+Fs],
      'the K2 curve reverse composite is the identity modulo the actual exceptional equation')
zero(Fs.subs({u:-s,v:1},simultaneous=True),'K2 has the regular whole-base z=0 section')
# All-rho S family:unit identities prove the parameterization is the WHOLE curve.
T=u*v+rho;Br=u*T*T-u+u*v+rho;E=Br-s
zero(E-(u*(T*T-1)+T-s),'the complete S_rho plane equation in its regular T parameter')
zero(s*s-1-(T*T-1)*(1+u*(s+T))+(s+T)*E,
     'T^2-1 is a unit on the whole plane curve over s!=+-1')
zero(Br-rho-u*(T*T+v-1),'u is a unit on the whole plane curve over s!=rho')
scurve={u:(s-t)/(t*t-1),v:(t-rho)*(t*t-1)/(s-t)}
zero(E.subs(scurve,simultaneous=True),'S_rho whole relative curve parameterization for arbitrary rho')
zero(T.subs(scurve,simultaneous=True)-t,'S_rho regular curve inverse t=uv+rho')
zero(((s-T)/(T*T-1)-u)*(T*T-1)+E,
     'S_rho first reverse composite holds in the whole plane-curve coordinate ring')
zero(((T-rho)*(T*T-1)/(s-T)-v)*(s-T)-v*E,
     'S_rho second reverse composite holds in the whole plane-curve coordinate ring')
zero(Br.subs({u:-s,v:rho/s},simultaneous=True)-s,'the t0=0 section is regular on S and works for every rho')
zero(T.subs({u:-s,v:rho/s},simultaneous=True),'the regular S section is exactly the pointed curve parameter zero')
for endpoint in [-1,1]:
    check((t-endpoint).subs(t,0)!=0,'the chosen t0 avoids the fixed puncture '+str(endpoint))
    zero(E.subs(T,endpoint)-(endpoint-s),'a fixed puncture has no additional points when its excluded value is inverted '+str(endpoint))
check((t-s).subs(t,0)==-s,'the sole extra s=0 deletion makes t0 avoid the moving puncture')
# General hyperbolic suspension:the x-open chart and zero fibre are literal.
pp=sp.Function('p')(z,w);susp=x*y-pp
zero(susp.subs(y,pp/x),'the generic x chart is the trivial A2 bundle over x!=0')
zero(susp.subs(x,0)+pp,'the whole special fibre is the affine y-line times p=0')
check(sp.diff(susp,y)==x,'smoothness away from x=0 is immediate from the unit y derivative')
zero(susp.subs({y:0,pp:0},simultaneous=True),'a section of p=0 supplies the entire x-line section')
algebra=len(checks);context=json.loads((HERE/'sources/motivic-family-barrier-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected full primary bytes '+entry['retained_file'])
check(len(context['newly_imported_proved_external_statements'])==3
      and all(e['hypotheses_checked'] for e in context['newly_imported_proved_external_statements']),
      'three explicit proved homotopy imports have checked smoothness,section and splitting hypotheses')
check(context['newly_used_unproved_external_statements']==[],'no conjectural homotopy or scheme-isomorphism classification is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==5 and not flags['native_smooth_fibre_has_required_additive_action']
      and not flags['native_family_identified_with_cotangent_family_as_schemes']
      and not flags['smooth_value_torsor_class_trivialized'] and not flags['critical_value_gluing_classified']
      and not flags['original_goal_complete'],
      'relative homotopy equivalence is not upgraded to an action,torsor trivialization,scheme isomorphism or full-source solution')
proof=HERE/'SMOOTH_VALUE_FAMILIES_AND_THE_MOTIVIC_HOMOTOPY_BARRIER.md'
source_paths=['sources/motivic-family-barrier-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in source_paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-motivic-family-barrier.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
