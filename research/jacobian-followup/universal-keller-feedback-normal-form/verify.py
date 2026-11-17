#!/usr/bin/env python3
"""Selected exact controls for the feedback normal form and formal source gauges."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,v):
    assert v,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def main():
    begun=time.monotonic();a,X=s.symbols('a X');v4=s.Matrix([-2*a*a*X**2,-2*X**4]);v5=s.Matrix([-s.Rational(5,2)*a*a*X**3,-3*X**5]);M=s.Matrix.hstack(v4,v5)
    ck('native basis wedge exact',M.det()==a*a*X**7)
    ck('native basis adjugate gives polynomial regularity test',M.adjugate()*M==a*a*X**7*s.eye(2))
    for j in [6,7,9]:
        vv=s.Matrix([-s.Rational(j,2)*a*a*X**(j-2),-(j-2)*X**j]);ck('native index '+str(j)+' exact two-control compression',all(s.expand(q)==0 for q in vv-(5-j)*X**(j-4)*v4-(j-4)*X**(j-5)*v5))
    c4,c5,c6,c7=s.symbols('c4 c5 c6 c7');coeff={4:c4,5:c5,6:c6,7:c7};AA=sum((5-j)*v*X**(j-4) for j,v in coeff.items());BB=sum((j-4)*v*X**(j-5) for j,v in coeff.items() if j>=5);raw=sum((s.Matrix([-s.Rational(j,2)*a*a*X**(j-2),-(j-2)*X**j])*v for j,v in coeff.items()),s.zeros(2,1))
    ck('arbitrary selected coefficient list compresses literally',all(s.expand(q)==0 for q in raw-M*s.Matrix([AA,BB])))
    ck('degree-six constant profile has unique quartic/quintic controls',AA.subs({c4:0,c5:0,c6:1,c7:0})==-X**2 and BB.subs({c4:0,c5:0,c6:1,c7:0})==2*X)
    lam=s.symbols('lambda');denom=1+s.Rational(4,3)*lam*a*a*X**6
    ck('formal pole is coprime to native basis degeneracy',s.expand(denom-a*X*(s.Rational(4,3)*lam*a*X**5))==1)
    # Exact full two-control determinant equation, with arbitrary symbolic controls.
    S,r=s.symbols('S r');A=s.Function('A')(S,r);B=s.Function('B')(S,r);V=s.Function('V')(S);W=s.Function('W')(S);h=s.Function('h')(S)
    Bout=(r-s.diff(h,S)-A*s.diff(V,S)-B*s.diff(W,S))/(2*S);Cout=r*S+2*h-S*s.diff(h,S)+A*(2*V-S*s.diff(V,S))+B*(2*W-S*s.diff(W,S))
    TT=S*s.diff(h,S,2)-s.diff(h,S)+A*(S*s.diff(V,S,2)-s.diff(V,S))+B*(S*s.diff(W,S,2)-s.diff(W,S));om=V*s.diff(W,S)-s.diff(V,S)*W
    rhs=V*s.diff(A,S)+((V/S-s.diff(V,S))*r+V*TT/S)*s.diff(A,r)+W*s.diff(B,S)+((W/S-s.diff(W,S))*r+W*TT/S)*s.diff(B,r)-om*(s.diff(A,S)*s.diff(B,r)-s.diff(A,r)*s.diff(B,S))
    actual=s.Matrix([Bout,Cout]).jacobian([S,r]).det()/(-r/S)
    ck('generic two-control nonlinear determinant exact',zero(r*(actual-1)-rhs))
    ck('generic response cross term equals native basis wedge',zero((om/S).subs({V:a**4*S**4,W:a**5*S**5}).doit()-a**9*S**7))
    # Laurent Wronskians and the source-jet exclusion, not a bounded classification.
    for n in [1,2,3]:
        pp=S**(-n)/s.Integer(2*n);qq=S**n;ck('Laurent power Wronskian '+str(n),zero(s.diff(pp,S)*qq-pp*s.diff(qq,S)+1/S))
    y,kappa,db2,dc4=s.symbols('y kappa db2 dc4');p0=1/(2*S);q0=S;b0=-1/(2*S)-s.Rational(3,2)*a*S;c0=S-a*S**3
    pp=p0+kappa*S;bb=b0-kappa*S+db2*S*S;cc=c0+dc4*S**4
    ck('native jets force remaining affine-r coefficient to zero',s.expand(s.diff(bb,S)*q0-pp*s.diff(cc,S)).coeff(S,1)==-2*kappa)
    ck('higher Laurent power has a genuine entire-source pole',zero((S**-2/s.Integer(4))*(-S*y)+y/(4*S)))
    f2,f3,f5=s.symbols('f2 f3 f5');FF=f2*S*S+f3*S**3+f5*S**5;HH=f2*S**4+s.Rational(6,5)*f3*S**5+s.Rational(10,7)*f5*S**7
    ck('retained-profile self equation integrates exactly',zero(s.diff(HH,S)-2*S*S*s.diff(FF,S)))
    # Entire cubic and stable signal normalizer identities.
    x,y,z,eta=s.symbols('x y z eta');xyz=[x,y,z];constructor=ROOT/'jacobian-followup/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(constructor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure source constructor only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(constructor),'exec'),env);compiler=env['source_map'];K=s.Matrix(compiler(3,{}));JK=K.jacobian(xyz);t=1+x*y;q=t*t*z+y*y*(1+3*t);sig=t*q;XX=x*q
    ck('actual cubic determinant one',zero(JK.det()-1));inv=JK.adjugate()
    ck('signal gradient is polynomial unimodular',zero((s.Matrix([sig]).jacobian(xyz)*inv[:,0])[0]-1))
    ck('exceptional components are comaximal by explicit Bezout',s.expand(x*x*q+t*(-3*t*t+5*t-1-t*x*x*z))==1)
    scaled={x:eta*x,y:y/eta,z:z/eta**2};KS=K.subs(scaled,simultaneous=True)
    ck('signal normalizer weighted scaling exact',zero(KS[0]-sig/eta**2))
    ck('compiler source scaling is old-output recoverable',all(zero(v) for v in KS-s.diag(eta**-2,eta**-1,eta)*K))
    ck('weighted base Jacobian is eta inverse squared',s.det(s.Matrix([eta*x,y/eta,z/eta**2]).jacobian(xyz))==eta**-2)
    alpha=s.symbols('alpha');swap_res=alpha/x**2-x**2/eta**2
    ck('swapped source roots fail t-zero polynomial landing',s.Poly(s.expand(swap_res*x*x),x).coeff_monomial(x**4)==-eta**-2)
    ss,gamma=s.symbols('s gamma');tt=1+x*y
    ck('generic-leaf t critical differential vanishes at distinguished cylinder',[s.diff(tt,v).subs({x:0,y:0}) for v in [x,y]]==[0,0])
    ck('swapped exceptional factor would have nonconstant critical value',(gamma*ss/tt).subs({x:0,y:0})==gamma*ss)
    ck('nonzero-signal fibre z reconstruction is exact',zero(sig.subs(z,(ss-t*y*y*(1+3*t))/t**3)-ss))
    w1,w2=s.symbols('w1 w2');stable=s.Matrix([eta*x,y/eta,z/eta**2,w1+q*w2**2,w2]);Jstable=stable.jacobian(xyz+[w1,w2])
    ck('actual stable normalizer has the stated vertical determinant',zero(Jstable.det()-eta**-2))
    # Formal source coordinate: exact first-order lift and selected second-order points.
    R=XX**3*(q+sig/3)-s.Rational(3,5)*XX**5;DD=sig**2*XX**6;vv=s.Matrix([0,-2*sig**2*XX**2,-2*XX**4]);tangent=R*vv;W0=inv*vv;A1=R*W0
    ck('formal first-order source lift is exact globally',all(zero(v) for v in JK*A1-tangent))
    ck('formal first-order source lift preserves signal',zero((s.Matrix([sig]).jacobian(xyz)*A1)[0]))
    points=[[1,0,1],[0,0,-1],[-s.Rational(1,2),3,26]];hess=[s.hessian(f,xyz) for f in K]
    for ii,point in enumerate(points):
        subs=dict(zip(xyz,point));Jp=JK.subs(subs);aa=A1.subs(subs);target2=(-s.Rational(4,3)*R*DD*vv).subs(subs);quad=s.Matrix([(aa.T*HH.subs(subs)*aa)[0]/2 for HH in hess]);a2=Jp.inv()*(target2-quad)
        ck('actual second-order source lift Taylor equation point '+str(ii),all(zero(v) for v in Jp*a2+quad-target2))
        ck('actual first-order lift signal tangent point '+str(ii),(s.Matrix([sig]).jacobian(xyz).subs(subs)*aa)[0]==0)
    vzero={x:0,y:0,z:-1};vone={x:-s.Rational(1,2),y:3,z:26}
    ck('polynomial-coordinate obstruction uses same actual old outputs',list(K.subs(vzero))==list(K.subs(vone))==[-1,0,0])
    ck('actual scalar tangent is not output-recoverable',list(tangent.subs(vzero))==[0,0,0] and list(tangent.subs(vone))==[0,s.Rational(32,15),s.Rational(32,15)])
    I=lambda T:T*(3*T-1)*(4*T+1);ck('marked multiplier retains the same non-descent witness',I(s.Integer(-1))==-12)
    RF=a**4*S**3*(r+s.Rational(1,3)-s.Rational(3,5)*a*S*S)
    ck('scalar formal tangent has genuine effective r dependence',s.diff(RF,r)==a**4*S**3)
    ck('degree-six compressed tangent still has effective r dependence',s.diff(-(a*S)**2*RF,r)!=0 and s.diff(2*a*S*RF,r)!=0)
    ck('all source-lift controls use exact rational arithmetic',not any(v.has(s.Float) for v in tangent))
    # An actual global polynomial family can jump degree despite formal source triviality.
    hp=S+a*S**3+lam*I(a)**2*a**4*S**4
    ck('formal marked equivalence is compatible with actual generic degree jump',s.degree(hp,S)==4 and s.degree(hp.subs(lam,0),S)==3)
    actualprofile=s.Matrix(compiler(4,{4:lam*I(sig)**2}));point0={x:1,y:0,z:1}
    ck('actual marked degree-jump family has exact Jacobian at a generic source point',actualprofile.jacobian(xyz).subs(point0).det()==1)
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local source pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes();ck('primary source pin '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),'effective_native_indices':[4,5],'generic_degree_notice':'These native indices do not bound the generic covering degree or total coordinate degree. The degree-six example is literally unchanged.','assurance':'Written universal module identities, characteristic-zero affine-r and stable signal-field rigidity, and a classical formal-etale source-gauge consequence. Selected exact controls are not independent whole-proof replay, formal verification, external review or publication-priority evidence. No new canonical atom or general multi-block classification.','sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
