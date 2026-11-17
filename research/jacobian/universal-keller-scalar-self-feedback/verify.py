#!/usr/bin/env python3
"""Selected exact controls for scalar self-feedback and formal integration."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,v):
    assert v,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def main():
    begun=time.monotonic();S,r,a,lam=s.symbols('S r a lambda');G=s.Function('g')(S,r);h=s.Function('h')(S);V=s.Function('V')(S)
    B=(r-s.diff(h,S)-G*s.diff(V,S))/(2*S);C=r*S+2*h-S*s.diff(h,S)+G*(2*V-S*s.diff(V,S));A=S*s.diff(h,S,2)-s.diff(h,S);Z=S*s.diff(V,S,2)-s.diff(V,S);L=S*s.diff(V,S)/V
    E=S*s.diff(G,S)+((1-L)*r+A+Z*G)*s.diff(G,r)
    ck('generic nonlinear response determinant exact',zero(s.Matrix([B,C]).jacobian([S,r]).det()/(-r/S)-1-V/(r*S)*E))
    g0=s.symbols('g0');Bf=(r-s.diff(h,S)-g0*s.diff(V,S))/(2*S);Cf=r*S+2*h-S*s.diff(h,S)+g0*(2*V-S*s.diff(V,S))
    ck('generic frozen source-chart determinant exact',zero(s.Matrix([Bf,Cf]).jacobian([S,r]).det()+r/S))
    # Reciprocal Riccati and degree/divisibility bookkeeping.
    vv=a**4*S**4+a**5*S**5;BB=S*s.diff(vv,S,2)-s.diff(vv,S);WW=s.integrate(s.expand(vv*BB/S**2),S)
    ck('mixed direction primitive derivative identity',zero(s.diff(WW,S)-vv*BB/S**2))
    ck('mixed direction lowest primitive order is six',min(m[0] for m,cf in s.Poly(WW,S).terms())==6)
    ck('mixed direction primitive degree is eight',s.degree(WW,S)==8)
    PP=s.cancel(vv/S**4);QQ=s.cancel(WW/S**6)
    ck('mixed direction normalized denominator degree exceeds numerator',s.degree(QQ,S)==2 and s.degree(PP,S)==1)
    C0=s.symbols('C0');ar=vv/(S*(C0+WW));Ric=S*s.diff(ar,S)+(1-S*s.diff(vv,S)/vv)*ar+BB*ar**2
    ck('general reciprocal Riccati family exact in mixed direction',zero(Ric))
    # Complete monomial exceptional branch, including a nontrivial background.
    bg,bzero=s.symbols('background bzero');j=4;v=a**j;vm=v*S**j;h0=S+a*S**3+bg*a**6*S**6;AA=S*s.diff(h0,S,2)-s.diff(h0,S);ZZ=S*s.diff(vm,S,2)-s.diff(vm,S);LL=j
    cc=s.Rational(j*(j-2),2*(j-1));a0=1/(cc*v);bpart=-a0*sum(cf*S**(n[0]+1-j)/s.Integer(n[0]+j-1) for n,cf in s.Poly(AA,S).terms());exception=a0*S**(1-j)*r+bpart+bzero*S**(2-2*j)
    op=lambda f:s.cancel(S*s.diff(f,S)+((1-LL)*r+AA+ZZ*f)*s.diff(f,r))
    ck('monomial exceptional Laurent branch solves full nonlinear equation',zero(op(exception)))
    ck('exceptional homogeneous term has strictly worse source pole',2-2*j<1-j)
    nohom=exception.subs(bzero,0).subs(r,1)
    ck('actual y-zero plane exceptional pole coefficient nonzero',zero(s.limit(S**(j-1)*nohom,S,0)-a0*j/s.Integer(j-1)) and a0!=0)
    ck('background derivative has required constant term',AA.subs(S,0)==-1)
    # A nonunit rational-control determinant would force the degenerate value zero.
    ee,vlead,alpha=s.symbols('e vlead alpha',nonzero=True)
    leading=-2*(ee-1)*alpha+ee*(ee-2)*vlead*alpha**2
    ck('whole-map leading regularity nonzero branch forces zero determinant',zero((1+vlead*leading).subs(alpha,1/(ee*vlead))))
    ck('exceptional rational control still has a genuine B-output pole',zero(-s.Rational(j,2)*v*(a0*j/s.Integer(j-1))+s.Rational(j,j-2)))
    # Entire source formulas and unperturbed-first-integral false lead.
    x,y,z=s.symbols('x y z');xyz=[x,y,z];constructor=ROOT/'jacobian/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(constructor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure compiler constructor only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(constructor),'exec'),env);compiler=env['source_map'];K=s.Matrix(compiler(3,{}));t=1+x*y;q=t*t*z+y*y*(1+3*t);sig=t*q;X=x*q
    ck('actual cubic Jacobian one',zero(K.jacobian(xyz).det()-1))
    reconx=S/r;recony=(1-r)/S;reconz=a*r**3-r*(r+3)*(1-r)**2/S**2
    ck('entire-source reconstruction has only Laurent chart denominators',zero(t.subs({x:reconx,y:recony},simultaneous=True)-1/r))
    ck('source reconstructed signal exact',zero(sig.subs({x:reconx,y:recony,z:reconz},simultaneous=True)-a))
    R=X**(j-1)*(q+sig/s.Integer(j-1))-s.Rational(3,j+1)*X**(j+1);D=sig**2*X**(2*j-2)
    RF=a**j*S**(j-1)*(r+s.Rational(1,j-1)-s.Rational(3,j+1)*a*S*S);DF=a**(2*j)*S**(2*j-2)
    ck('polynomial first integral reconstruction exact',zero(R.subs({x:reconx,y:recony,z:reconz},simultaneous=True)-RF))
    ck('polynomial denominator profile reconstruction exact',zero(D.subs({x:reconx,y:recony,z:reconz},simultaneous=True)-DF))
    baseA=3*a*S*S-1;respB=j*(j-2)*a**j*S**(j-1)
    old=lambda f:s.cancel(S*s.diff(f,S)+((1-j)*r+baseA)*s.diff(f,r))
    nonlinear=lambda f:s.cancel(S*s.diff(f,S)+((1-j)*r+baseA+respB*f)*s.diff(f,r))
    ck('naive R is an exact first integral only of frozen response',zero(old(RF)))
    defect=vm/(r*S)*nonlinear(lam*RF);expected=lam**2*j*(j-2)*RF*a**3*(a*S)**(3*j-3)/r
    ck('naive self-feedback has exact quadratic determinant defect',zero(defect-expected) and defect!=0)
    point={x:s.Integer(1),y:s.Integer(0),z:s.Integer(1)};literal=s.Matrix(compiler(j,{j:lam*R}));Jpoint=literal.jacobian(xyz).subs(point)
    ck('actual naive whole-source Jacobian has quadratic defect',zero(Jpoint.det()-1-lam**2*s.Rational(88,15)))
    ck('actual naive scalar control is not a retained signal',R.subs({y:0,z:1}).has(x))
    for name,gval in [('t',t),('q',q)]:
        actual=s.Matrix(compiler(j,{j:gval}));detpoint=actual.jacobian(xyz).subs(point).det();ck('actual additional source channel '+name+' is not Keller',detpoint!=1)
    u=s.symbols('u');positive=s.Matrix(compiler(j,{j:sig+u}))
    ck('actual retained-signal control Jacobian point one',positive.jacobian(xyz).subs(point).det()==1)
    ck('all parameter response components vanish at x zero',all(f.subs(x,0)==0 for f in [-s.Rational(j,2)*t*t*x**(j-2)*q**j,-(j-2)*X**j]))
    ck('rational retained control may cancel a common direction factor',s.cancel(sig*(1/sig))==1)
    # Exact rational solution and every finite-order polynomial truncation.
    gformal=lam*RF/(1+cc*lam*DF)
    ck('rational formal sum solves exact nonlinear Keller equation',zero(nonlinear(gformal)))
    Ecurrent=lambda f,g:s.cancel(S*s.diff(f,S)+((1-j)*r+baseA+respB*g)*s.diff(f,r))
    trial=s.Function('trial')(S,r);constraint=(1+cc*lam*DF)*trial-lam*RF
    ck('universal implicit constraint derivative identity',zero(Ecurrent(constraint,trial)-(1+cc*lam*DF)*Ecurrent(trial,trial)))
    for n in [1,2,3]:
        trunc=sum((-cc)**k*lam**(k+1)*RF*DF**k for k in range(n));res=s.Poly(s.cancel(nonlinear(trunc)).expand(),lam)
        ck('parameter truncation '+str(n)+' solves modulo next power',all(m[0]>=n+1 for m,cf in res.terms() if cf!=0))
        ck('parameter truncation '+str(n)+' geometric constraint exact',zero((1+cc*lam*DF)*trunc-lam*RF+lam*RF*(-cc*lam*DF)**n))
    g2=lam*R-cc*lam**2*R*D;actual2=s.Matrix(compiler(j,{j:g2}));det2=s.factor(actual2.jacobian(xyz).subs(point).det()-1)
    ck('literal second-order polynomial map is Keller mod lambda cubed',s.rem(det2,lam**3,lam)==0)
    ck('literal second-order polynomial map is not exactly Keller',det2!=0)
    ck('literal second-order coefficients are exact rationals',not any(f.has(s.Float) for f in actual2))
    # Genuine whole-source pole, detected already on y=0,z=1.
    denom=1+cc*lam*D;pole={x:s.Integer(1),y:s.Integer(0),z:s.Integer(1),lam:-1/cc}
    ck('rational-sum denominator has actual source zero',denom.subs(pole)==0)
    ck('rational-sum last-output numerator does not cancel that pole',(-(j-2)*X**j*lam*R).subs(pole)!=0)
    sliceQ=denom.subs({y:0,z:1});sliceN=lam*R.subs({y:0,z:1})
    ck('rational source pole survives polynomial gcd check',s.gcd(sliceQ,sliceN)==1)
    # Retained marked planes/thickenings and canonical idempotent fibre relation.
    I=lambda T:T*(3*T-1)*(4*T+1);mark=I(a)**2;marked=lam*mark*RF/(1+cc*lam*mark*DF)
    ck('marked rational sum still solves exact Keller equation',zero(nonlinear(marked)))
    for k in [0,1,2]:
        coeff=(-cc)**k*mark**(k+1)*RF*DF**k;ck('formal marked coefficient '+str(k)+' divisible by full thickening',zero(coeff/mark-(-cc)**k*mark**k*RF*DF**k))
    ck('three complete marked-plane first values annihilate I',[I(v) for v in [s.Integer(0),s.Rational(1,3),-s.Rational(1,4)]]==[0,0,0])
    beta=s.symbols('beta');ck('canonical whole-fibre idempotent ideal identity',zero(I(beta/3)-beta*(beta-1)*(4*beta+3)/9))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local source pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes();ck('primary source pin '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),'finite_parameter_orders_controlled':[1,2,3],'scope':'Written all-degree Laurent classification, whole-source scalar-line rigidity and formal marked integration boundary. The finite-order examples are selected exact controls, not an inference of the universal statements. The quasi-translation and PDE sources are contextual, not imported classification dependencies.','sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Not independent whole-proof replay, formal verification, external review or publication-priority evidence. No new canonical atomic map, degree bound or whole multi-block classification.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
