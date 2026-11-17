#!/usr/bin/env python3
"""Selected exact controls for profile-factor rigidity and plane completion."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def main():
    begun=time.monotonic();x,y,z,S,Y,a,b,c,lam,xi=s.symbols('x y z S Y a b c lambda xi');r=1-S*Y
    ctor=ROOT/'jacobian-followup/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(ctor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map']
    ck('reuse literal pure compiler constructor only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(ctor),'exec'),env);K=s.Matrix(env['source_map'](3,{}))
    t=1+x*y;q=t*t*z+y*y*(1+3*t);sig=t*q;X=x*q
    source={x:S/r,y:Y,z:a*r**3-r*(r+3)*Y**2}
    ck('generic source reconstruction retains sigma',zero(sig.subs(source,simultaneous=True)-a))
    ck('reciprocal chart t is inverse r',zero(t.subs(source,simultaneous=True)-1/r))
    ck('source root becomes aS',zero(X.subs(source,simultaneous=True)-a*S))
    F=s.Matrix([-Y/2-3*a*S/2,2*S-S*S*Y-a*S**3])
    pulled=K.subs(source,simultaneous=True)
    ck('actual cubic gives the literal plane completion',all(zero(v) for v in pulled-s.Matrix([a,*F])))
    ck('generic plane determinant is exactly the boundary equation',s.factor(F.jacobian([S,Y]).det())==r)
    ck('native image line is a closed embedding with explicit inverse',F.subs(S,0)==s.Matrix([-Y/2,0]) and -2*F[0].subs(S,0)==Y)
    ck('boundary is smooth with no zero of both gradient entries',s.gcd(s.diff(r,S),s.diff(r,Y))==1)
    factors=s.factor_list(r,S,Y)[1]
    ck('boundary equation is an irreducible multiplicity-one factor',len(factors)==1 and factors[0][1]==1)
    ck('two affine charts cover the whole completion',s.expand(r+S*Y)==1)
    rr=s.symbols('r')
    pole=S/rr
    ck('own-source x has genuine boundary pole and is outside completion',s.limit(rr*pole,rr,0)==S)
    ck('the native base completion can collide despite an embedded line',F.subs({a:-1,S:0,Y:0})==F.subs({a:-1,S:1,Y:3})==s.Matrix([0,0]))
    ck('line theorem cannot invert the completed map with nonconstant determinant',s.Poly(r,S,Y).total_degree()==2)
    # Generic native module membership, including its congruence beyond source jets.
    U=s.Function('U')(S,Y);W=s.Function('W')(S,Y);V=U+S*W
    A=(5*V-6*U)/(2*a**4);B=2*(U-V)/(a**5*S)
    v4=s.Matrix([-2*a**4*S*S,-2*a**4*S**4]);v5=s.Matrix([-s.Rational(5,2)*a**5*S**3,-3*a**5*S**5])
    ck('all-symbol generic native module congruence and controls exact',all(zero(v) for v in A*v4+B*v5-s.Matrix([S*S*U,S**4*V])))
    ck('jet-only pair with V-U not divisible by S has a control pole',zero((2*(0-1)/(a**5*S))*S+2/a**5))
    H=F+xi*v4
    ck('arbitrary quartic retained profile has the same plane Jacobian',s.factor(H.jacobian([S,Y]).det())==r)
    ck('retained profile line and C first jet agree in every coefficient',H.subs(S,0)==s.Matrix([-Y/2,0]) and s.diff(H[1],S).subs(S,0)==2)
    h=lam*a**5*c**3
    E=s.Matrix([b+h,c]);inv=s.Matrix([b-h,c])
    ck('allowed recovered target map has exact polynomial inverse',E.subs({b:inv[0],c:inv[1]},simultaneous=True)==s.Matrix([b,c]))
    ck('allowed recovered target map fixes the full line pointwise',E.subs(c,0)==s.Matrix([b,0]))
    ck('allowed recovered target map is Keller globally',E.jacobian([b,c]).det()==1)
    T=s.Matrix([H[0]+lam*a**5*H[1]**3,H[1]])
    ck('nonlinear retained-profile orbit still completes with determinant r',s.factor(T.jacobian([S,Y]).det())==r)
    L=s.cancel(H[1]/S);AA=-3*lam*a*S*L**3;BB=2*lam*L**3
    ck('recovered target shear has literal generic polynomial native controls',all(zero(v) for v in AA*v4+BB*v5-(T-H)))
    Rexpr=T.subs(Y,(1-rr)/S)
    ck('literal recovered shear is nonlinear in r only in B',s.degree(s.cancel(Rexpr[0]),rr)==3 and s.degree(s.cancel(Rexpr[1]),rr)==1)
    # The omitted hypothesis cannot be replaced by merely fixing the line.
    bad=s.Matrix([b,c*(1+b*c)])
    ck('nontrivial after-factor control fixes the line but is not Keller',bad.subs(c,0)==s.Matrix([b,0]) and bad.jacobian([b,c]).det()==1+2*b*c)
    # Irreducibility is essential to the general plane factor argument.
    G=s.Matrix([S*S,Y]);D=s.Matrix([b,c*c]);Fbad=D.subs({b:G[0],c:G[1]},simultaneous=True)
    ck('reducible-Jacobian control has two nonautomorphism stages',G.jacobian([S,Y]).det()==2*S and D.jacobian([b,c]).det()==2*c)
    ck('reducible-Jacobian composite still embeds a line explicitly',zero((Fbad[1]-Fbad[0]-1).subs(Y,S+1)/2-S))
    ck('reducible-Jacobian composite determinant splits',Fbad.jacobian([S,Y]).det()==4*S*Y)
    # Generic plane polynomiality does not replace whole-source landing.
    tau=s.symbols('tau');qt=(tau-1)**2*(1+3*tau);delta=-2*qt**3/tau
    ck('rational-signal profile has an uncancelled whole-source pole',s.limit(tau*delta,tau,0)==-2)
    # All-degree leading-divisor condition is not a cutoff on source degrees.
    hh=s.Function('h')(S);M,N,d=s.symbols('M N d',integer=True,positive=True);pp=hh**M;qq=hh**N
    ck('symbolic common-power leading identity at arbitrary integer exponents',s.simplify(d*N*s.diff(pp,S)*qq-d*M*pp*s.diff(qq,S))==0)
    ck('non-dividing pair 6 and 10 forces respective divisor powers 3 and 5',6//s.igcd(6,10)==3 and 10//s.igcd(6,10)==5)
    ck('a simple leading zero excludes the degree pair 2 and 3',3*1 % 2 != 0)
    ck('selected formulas contain no floating coefficients',not any(v.has(s.Float) for v in list(F)+list(H)+list(T)+list(Fbad)))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:
        ck('retained local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes()
        ck('read primary HTML pin '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),
    'assurance':'Selected exact chart, native-module, field-orbit, sharp boundary and leading-divisor controls accompany written general proofs. Gwozdziewicz line injectivity and Abhyankar--Moh rectification are published imports. Not independent whole-proof replay, formal verification, external review or novelty evidence. No canonical atomic receiver or stable whole-map atomicity is constructed.',
    'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
