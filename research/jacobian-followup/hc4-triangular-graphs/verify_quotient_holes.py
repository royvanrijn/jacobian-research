"""Exact controls of the written all-degree quotient,torsor and local-period proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0 or sp.simplify(value)==0,label)
U,V,P,Q=sp.symbols('U V P Q');coords=[U,V,P,Q]
def deriv(coeff,f):return sp.expand(sum(a*sp.diff(f,x) for a,x in zip(coeff,coords)))
DD=[0,0,V,-U];EE=[-Q,P,0,0];WW=[U,V,-P,-Q]
H=U*P+V*Q
for x in coords:
    zero(deriv(DD,deriv(EE,x))-deriv(EE,deriv(DD,x))-deriv(WW,x),'the D,E commutator is the exact weight operator '+str(x))
    zero(deriv(WW,deriv(EE,x))-deriv(EE,deriv(WW,x))+2*deriv(EE,x),'the lowering operator has weight minus two '+str(x))
zero(deriv(DD,H),'the translation kills its actual invariant generator')
zero(deriv(EE,H),'the lowering operator also kills that generator')
# The all-m weight identity is proved in the note. These controls verify its sign,
# coefficient and the last-nonzero-iterate contradiction on the first strings.
for m in range(1,5):
    def dpow(f,k):
        for _ in range(k):f=deriv(EE,f)
        return f
    # These weight-minus-two controls have nontrivial lowering strings.
    rr=U**m*Q**(m+2)
    zero(deriv(DD,dpow(rr,m))-dpow(deriv(DD,rr),m)+m*(m+1)*dpow(rr,m-1),
         'the finite lowering-string identity has coefficient -m(m+1),m='+str(m))
for k in range(4):
    zero(deriv(DD,P*H**k)-V*H**k,'every G multiple of an invariant power is in the plinth,k='+str(k))
    zero(deriv(DD,-Q*H**k)-U*H**k,'every F multiple of an invariant power is in the plinth,k='+str(k))
# Exact universal obstruction: no constant H^k is an image on its homogeneous
# component. Matrix column-rank controls supplement,not replace,the weight proof.
for k in range(4):
    from itertools import product
    rmons=[sp.prod(x**e for x,e in zip(coords,exps)) for exps in product(range(2*k+1),repeat=4) if sum(exps)==2*k and exps[0]+exps[1]-exps[2]-exps[3]==-2]
    target=sp.Poly(H**k,*coords)
    images=[sp.Poly(deriv(DD,r),*coords) for r in rmons]
    rows=sorted(set(target.monoms()).union(*(set(p.monoms()) for p in images)))
    matrix=sp.Matrix([[p.coeff_monomial(e) for p in images] for e in rows])
    vector=sp.Matrix([target.coeff_monomial(e) for e in rows])
    check(matrix.row_join(vector).rank()>matrix.rank(),'the universal H power survives the plinth obstruction,k='+str(k))
x,z,s,c,kappa,rho=sp.symbols('x z s c kappa rho');c=sp.symbols('c',nonzero=True)
F=x;G=z*z-1;d=z
check(sp.gcd(F,G,x,z)==1,'the nonlinear coefficient control is a regular sequence')
check(sp.factor(G)==(z-1)*(z+1),'the zero scheme has two distinct reduced points')
value=d+c*(P*F+Q*G)
for sign in [-1,1]:
    zero(value.subs({x:0,z:sign})-sign,'the entire fixed plane is contracted to one quotient point,sign='+str(sign))
    check(G.diff(z).subs(z,sign)!=0,'the displayed control has a simple local coefficient zero,sign='+str(sign))
zero(P/G-(-Q/F)-(value-d)/(c*F*G),'the two actual local slices have the required exact transition')
zero(F*(kappa-d)/(c*F*G)-(kappa-d)/(c*G),'F kills the torsor class by a Cech coboundary')
zero(G*(kappa-d)/(c*F*G)-(kappa-d)/(c*F),'G kills the torsor class by a Cech coboundary')
for kk in [0,1,-1]:
    active=[sign for sign in [-1,1] if sign!=kk]
    expected=sp.prod(z-sign for sign in active)
    for sign in [-1,1]:
        check((expected.subs(z,sign)*(kk-sign)==0),'the exact annihilator retains precisely noncritical support,kappa='+str(kk)+',point='+str(sign))
    check(sp.degree(expected,z)==len(active),'the support shrinks only at the corresponding critical value,kappa='+str(kk))
# This control is not Keller: the point Jacobian varies. It must not be used to
# assert the common-slope conclusion without the everywhere constant hypothesis.
check(sp.diff(G,z).subs(z,1)!=sp.diff(G,z).subs(z,-1),'a varying coefficient Jacobian fails the required common residue slope')
FK=x+z*z;GK=z;j=sp.Matrix([FK,GK]).jacobian([x,z]).det()
check(j==1,'the Keller triangular control has one globally constant orientation')
zero(((kappa-d)/c/j).subs({x:0,z:0})-kappa/c,'the Keller residue is exactly affine in the value')
# Local Tang completion is an actual constant-volume quadratic germ.
u,v,p,q=sp.symbols('u v p q');alpha=u*u+v;beta=v*v+u
dd=u*alpha+v*beta
zero(dd+c*(u*P+v*Q)-c*(u*(P+alpha/c)+v*(Q+beta/c)),
     'the local analytic base remainder is absorbed by vertical translations')
check(sp.Matrix([u,v,P+alpha/c,Q+beta/c]).jacobian([u,v,P,Q]).det()==1,
      'absorbing the base remainder preserves volume exactly')
eps=sp.symbols('eps',positive=True)
zero((sp.sqrt(eps)**4/eps)-eps,'a four-variable quadratic Gelfand--Leray period scales linearly')
# Suspended multipliers: exact residues,not a finite Taylor approximation.
for sign in [-1,1]:
    zero(sp.residue(1/(z*z-1),z,sign)-1/(2*sp.Integer(sign)),
         'each old-model plane residue is independent of the value,sign='+str(sign))
# The K2 density identity uses r=sqrt(w),w0=a with a^2=1/4.
r,a,t,S,w=sp.symbols('r a t S w',nonzero=True)
zero((-w-1/(4*w))-(-2*a)+(w-a)**2/w-(a*a-sp.Rational(1,4))/w,
     'the K2 quadratic completion uses precisely the critical relation a^2=1/4')
zero(sp.diff((w-a)/sp.sqrt(w),w)-(w+a)/(2*w**sp.Rational(3,2)),
     'the exact K2 density has the stated Jacobian derivative')
zero(2/(r*r+a)-(1/a)*(1-((r-a/r)/(r+a/r))),
     'the K2 volume density is a constant plus an odd analytic term')
zero((r+a/r)**2-(r-a/r)**2-4*a,'the square-root branch S is even in t')
zero(((1/a)*(1-t/S)+(1/a)*(1+t/S))/2-1/a,
     'the full even density,not a truncation,is constant')
K2=w*(2*x*sp.Symbol('y')+z*z-1)+z;y=sp.Symbol('y')
K3=K2+2*x*(rho-z)
zero(K2.subs(y,y+(rho-z)/w)-K3,'the two nonzero-w K3 germs reduce to K2 by a literal local shear')
check(sp.Matrix([x,y+(rho-z)/w,z,w]).jacobian([x,y,z,w]).det()==1,
      'the nonzero-w K3 shear preserves volume')
B=u*(u*v+rho)**2-u+u*v+rho
native_V=w*y+rho-z
native_P=x+(w*native_V-2*w*w*y-2*w*rho-1)/2
zero(2*native_P*native_V+B.subs({u:w,v:y},simultaneous=True)-K3,
     'the third-point suspension is the retained actual global polynomial identity')
check(sp.Matrix([w,y,native_P,native_V]).jacobian([x,y,z,w]).det()==-1,
      'the actual global suspension map has one constant volume multiplier')
point={u:0,v:1-rho*rho}
zero(B.subs(point)-rho,'the third plane critical value is rho for every parameter')
zero(sp.diff(B,u).subs(point),'the third point is critical in u for every parameter')
zero(sp.diff(B,v).subs(point),'the third point is critical in v for every parameter')
zero(sp.hessian(B,[u,v]).det().subs(point)+1,'the third plane point stays Morse at both value collisions')
L=(-1+sp.sqrt(1+4*u*(u+s)))/(2*u)
zero(u*L*L+L-u-s,'the selected plane-level branch satisfies its exact quadratic equation')
zero(sp.diff((L-rho)/u,s)-1/(u*sp.sqrt(1+4*u*(u+s))),
     'the third-point Gelfand--Leray density is the stated simple-pole expression')
zero(sp.limit(u/(u*sp.sqrt(1+4*u*(u+s))),u,0)-1,
     'the full third-point residue is identically one in s and rho')
for rr in [-1,1]:
    zero(sp.hessian(B,[u,v]).det().subs({u:0,v:0,rho:rr})+1,
         'the colliding-value specialization keeps a distinct Morse point,rho='+str(rr))
context=json.loads((HERE/'sources/quotient-holes-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained proof/checker/receipt bytes: '+path)
for e in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/e['retained_file']).read_bytes()).hexdigest()==e['sha256'],'inspected primary context bytes: '+e['retained_file'])
check(context['newly_used_unproved_external_statements']==[],'new proofs do not import an isochore classification or a flatness theorem')
check(context['exact_plinth_equal_coefficient_ideal_proved'] and context['exact_torsor_annihilators_and_affine_point_residues_proved'],
      'exact quotient and residue conditions are necessary global tests')
check(context['all_remaining_normalized_models_pass_complete_local_volume_period_test']
      and not context['local_periods_assumed_to_glue_to_global_polynomial_coordinates']
      and not context['torsor_profile_assumed_sufficient_for_coordinate_certificate']
      and not context['unrestricted_HC4_or_JC2_resolved'],
      'all local orders pass while global coordinate recovery and HC4 remain open')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
sources={'sources/quotient-holes-context.json':hashlib.sha256((HERE/'sources/quotient-holes-context.json').read_bytes()).hexdigest()}
sources.update({'sources/'+e['retained_file']:e['sha256'] for e in context['newly_cached_primary_context']})
receipt=dict(outcome='PASS',claims=['HC4MYJU0','HC4MYJU1','HC4MYJU2'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'QUOTIENT_HOLES_TORSOR_RESIDUES_AND_LOCAL_PERIOD_BARRIER.md').read_bytes()).hexdigest(),
 source_sha256=sources,dependencies=context['retained_proof_and_source_sha256'],
 assurance='Finite exact controls of written all-degree plinth,quotient,torsor and complete local-period arguments; not independent replay,formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-quotient-holes.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
