#!/usr/bin/env python3
"""Selected exact embedded-curve and arbitrary-intermediate surface controls."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def main():
    begun=time.monotonic();x,y,z,S,Y,a,b,c,w,v,h,j,k,rr,dd=s.symbols('x y z S Y a b c w v h j k r d')
    ctor=ROOT/'jacobian/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(ctor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map']
    ck('reuse pure compiler constructor only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(ctor),'exec'),env);K=s.Matrix(env['source_map'](3,{}))
    t=1+x*y;q=t*t*z+y*y*(1+3*t);sig=t*q;X=x*q
    witness={x:0,y:Y,z:a-4*Y**2}
    ck('actual generic source witness retains signal',sig.subs(witness,simultaneous=True)==a)
    ck('actual witness image is an entire closed embedded target line',K.subs(witness,simultaneous=True)==s.Matrix([a,-Y/2,0]))
    ck('curve parameter has exact polynomial recovery',-2*K[1].subs(witness,simultaneous=True)==Y)
    v4=s.Matrix([-2*sig**2*X**2,-2*X**4]);v5=s.Matrix([-s.Rational(5,2)*sig**2*X**3,-3*X**5])
    ck('both universal native controls vanish on the entire witness',v4.subs(witness,simultaneous=True)==s.zeros(2,1) and v5.subs(witness,simultaneous=True)==s.zeros(2,1))
    # Nonlinear intermediate embeddings are enough; a straight curve is unnecessary.
    A=s.Matrix([b+(c+b*b)**2,c+b*b]);Ai=s.Matrix([b-c*c,c-(b-c*c)**2])
    ck('nonlinear plane rectifier has exact inverse',Ai.subs({b:A[0],c:A[1]},simultaneous=True)==s.Matrix([b,c]))
    ck('nonlinear rectifier and its inverse are Keller',A.jacobian([b,c]).det()==1 and Ai.jacobian([b,c]).det()==1)
    gamma=A.subs({b:Y,c:0},simultaneous=True)
    ck('intermediate curve can be nonlinear in both coordinates',gamma==s.Matrix([Y+Y**4,Y**2]))
    ck('coordinate-ring sandwich has explicit curve inverse',s.expand(gamma[0]-gamma[1]**2)==Y)
    ck('Keller target of nonlinear curve has embedded final line',Ai.subs({b:gamma[0],c:gamma[1]},simultaneous=True)==s.Matrix([Y,0]))
    # The restriction can be polynomial while the surrounding generic map has a pole.
    r=1-S*Y;M=s.Matrix([S/r,Y*r]);Mi=s.Matrix([S*r,Y/r])
    ck('generic surface automorphism retains boundary unit',zero(1-M[0]*M[1]-r))
    ck('generic pole-map has exact inverse on the open surface',all(zero(vv) for vv in Mi.subs({S:M[0],Y:M[1]},simultaneous=True)-s.Matrix([S,Y])))
    ck('generic pole-map is regular and embeds the whole carried line',M.subs(S,0)==s.Matrix([0,Y]))
    ck('the same map has a genuine pole elsewhere on the boundary',s.limit(rr*M[0].subs(Y,(1-rr)/S),rr,0)==S)
    # Arbitrary identity and signal mixing in an actual intermediate automorphism.
    mix=s.Matrix([w+b*b,a+c*c,b+(w+b*b)*(a+c*c),c])
    invmix=s.Matrix([v-(j-v*h)**2,h-k*k,j-v*h,k])
    ck('all-coordinate mixing has an exact inverse',invmix.subs({v:mix[0],h:mix[1],j:mix[2],k:mix[3]},simultaneous=True)==s.Matrix([w,a,b,c]))
    ck('all-coordinate mixing is globally Keller',mix.jacobian([w,a,b,c]).det()==1)
    ck('retained tuple fibre has an explicit whole plane parameterization',invmix[:2,0].subs({v:mix[0],h:mix[1],j:mix[2],k:mix[3]},simultaneous=True)==s.Matrix([w,a]))
    mixed_curve=mix.subs({b:-Y/2,c:0},simultaneous=True)
    ck('mixed intermediate carried curve still recovers its parameter',s.expand(-2*(mixed_curve[2]-mixed_curve[0]*mixed_curve[1]))==Y)
    ck('mixed intermediate curve keeps the generic retained values',invmix.subs({v:mixed_curve[0],h:mixed_curve[1],j:mixed_curve[2],k:mixed_curve[3]},simultaneous=True)==s.Matrix([w,a,-Y/2,0]))
    # A whole closed component does not remove the other inverse-line components.
    F=s.Matrix([-Y/2-3*a*S/2,2*S-S*S*Y-a*S**3]);D=2-S*Y-a*S*S
    ck('actual generic map is etale on the open source surface',s.factor(F.jacobian([S,Y]).det())==r)
    ck('target line preimage splits off the whole carried component',s.expand(F[1]-S*D)==0 and s.expand((D+S*(Y+a*S))/2)==1)
    otherY=(2-a*S*S)/S
    ck('remaining inverse-line component is distinct and explicit',zero(D.subs(Y,otherY)) and zero(F[0].subs(Y,otherY)+1/S+a*S))
    ck('other component is etale precisely off its deleted boundary points',zero(s.diff(-1/S-a*S,S)+(a*S*S-1)/S**2))
    ck('boundary unit is irreducible in the polynomial plane',len(s.factor_list(r,S,Y)[1])==1)
    # Etale open immersions alone do not obey whole-factor target saturation.
    rational=1/(1-b*c)
    ck('open-immersion control has a nonpolynomial rational target pullback',s.Poly(s.denom(s.cancel(rational)),b,c).total_degree()==2 and zero(rational.subs({b:S,c:Y},simultaneous=True)-1/r))
    # All-degree etale self-maps of the source surface; these are not whole lifts.
    d=s.symbols('d',integer=True,positive=True);HH=(1-rr**d)/(1-rr)
    ck('all-degree geometric-sum Jacobian identity',s.simplify(HH-(1-rr)*s.diff(HH,rr)-d*rr**(d-1))==0)
    ck('all-degree canonical volume ratio identity',s.simplify(rr**d*d*rr**(d-1)/rr-d*rr**(2*d-2))==0)
    for nd in [2,3]:
        Hn=sum(r**ell for ell in range(nd));Jn=s.Matrix([S,Y*Hn])
        ck('generic self-map unit power degree '+str(nd),zero(1-Jn[0]*Jn[1]-r**nd))
        ck('generic self-map etale Jacobian degree '+str(nd),s.factor(Jn.jacobian([S,Y]).det()-nd*r**(nd-1))==0)
        ck('generic degree is exact on the S-nonzero chart '+str(nd),s.degree(rr**nd-a,rr)==nd)
        ck('generic self-map retains the whole affine line '+str(nd),Jn.subs(S,0)==s.Matrix([0,nd*Y]))
        tau=s.symbols('tau');ylift=(tau-1)*sum(tau**(-ell) for ell in range(nd))
        ck('canonical whole-source lift has a genuine source pole degree '+str(nd),s.limit(tau**(nd-1)*ylift,tau,0)==-1)
    # Higher ambient dimension and a non-Keller target cannot replace plane cancellation.
    bad=s.Matrix([b,c*c]);ck('non-Keller target can keep a closed line and remain noninvertible',bad.subs({b:Y,c:0},simultaneous=True)==s.Matrix([Y,0]) and bad.jacobian([b,c]).det()==2*c)
    ck('three-dimensional boundary control is actually Keller',zero(K.jacobian([x,y,z]).det()-1))
    curve3=K.subs({x:0,y:Y,z:0},simultaneous=True)
    ck('three-dimensional noninvertible target has an embedded source curve',curve3==s.Matrix([4*Y**2,-Y/2,0]) and -2*curve3[1]==Y)
    ck('three-dimensional target is nevertheless noninjective',K.subs({x:0,y:0,z:-1})==K.subs({x:-s.Rational(1,2),y:3,z:26})==s.Matrix([-1,0,0]))
    ck('selected maps use exact coefficients only',not any(vv.has(s.Float) for vv in list(K)+list(mix)+list(F)))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:
        ck('retained local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes();ck('primary HTML pin '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),
    'assurance':'Selected exact witness, nonlinear curve, arbitrary intermediate coordinate mixing, open-surface pole and all-degree self-map identities accompany written general proofs. Published plane line injectivity, embedded-line rectification and etale base-change/section facts are imports. Not independent whole-proof replay, formal verification, external review or publication-priority evidence. No canonical atomic receiver or finite polynomial scalar-tangent integration is constructed.',
    'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
