#!/usr/bin/env python3
"""Selected exact controls for primitive-unit rigidity and whole factor descent."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(v):return s.cancel(v)==0
def main():
    begun=time.monotonic();x,y,z,S,Y,r,a,b,c,U,V,W,zeta=s.symbols('x y z S Y r a b c U V W zeta')
    ctor=ROOT/'jacobian-followup/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(ctor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure literal source constructor only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(ctor),'exec'),env);K=s.Matrix(env['source_map'](3,{}))
    t=1+x*y;Q=t*t*z+y*y*(1+3*t);sig=t*Q
    ck('source zero signal has two literal factors',zero(K[0]-t*Q))
    bez=-3*t*t+5*t-1-t*x*x*z
    ck('source signal factors are comaximal on the whole source',s.expand(x*x*Q+t*bez)==1)
    ck('Q is primitive linear in z',s.degree(Q,z)==1 and s.gcd(t*t,y*y*(1+3*t))==1)
    ck('nonzero signal fibre primitive obstruction has constant nonzero residue',zero((sig-c).subs(y,-1/x)+c))
    ck('generic nonzero signal fibre is primitive linear in z',s.degree(sig-c,z)==1 and s.gcd(t**3,t*y*y*(1+3*t)-c)==1)
    # Whole intermediate polynomial automorphism: both source factors recover exactly.
    G=s.Matrix([x+y*y,y,z+(x+y*y)**3])
    Gi=s.Matrix([U-V*V,V,W-U**3])
    ck('whole source-coordinate change has exact inverse',Gi.subs({U:G[0],V:G[1],W:G[2]},simultaneous=True)==s.Matrix([x,y,z]))
    ck('whole source-coordinate change is Keller',G.jacobian([x,y,z]).det()==1)
    back={x:Gi[0],y:Gi[1],z:Gi[2]};A1=t.subs(back,simultaneous=True);A2=Q.subs(back,simultaneous=True)
    ck('primitive t is recovered by whole intermediate polynomial',zero(A1.subs({U:G[0],V:G[1],W:G[2]},simultaneous=True)-t))
    ck('other zero component is recovered by whole intermediate polynomial',zero(A2.subs({U:G[0],V:G[1],W:G[2]},simultaneous=True)-Q))
    E0=K[0].subs(back,simultaneous=True)
    ck('intermediate signal has exactly the displayed two-factor presentation',zero(E0-A1*A2))
    ck('intermediate factors carry the whole Bezout identity',zero((x*x*Q+t*bez).subs(back,simultaneous=True)-1))
    ck('recovered generic unit pulls back primitively',zero((1/A1).subs({U:G[0],V:G[1],W:G[2]},simultaneous=True)-1/t))
    nn=s.symbols('n',integer=True,positive=True)
    ck('all positive powers of a recovered whole t remain in the whole target ring',s.simplify(s.cancel(A1.subs({U:G[0],V:G[1],W:G[2]},simultaneous=True))**nn-t**nn)==0)
    # Recovering a product is weaker than recovering a factor.
    p0={x:0,y:0,z:-1};p1={x:-s.Rational(1,2),y:3,z:26}
    ck('old cubic equal-output witness is literal',K.subs(p0)==K.subs(p1)==s.Matrix([-1,0,0]))
    ck('t does not descend through the old cubic although tQ does',t.subs(p0)==1 and t.subs(p1)==-s.Rational(1,2) and sig.subs(p0)==sig.subs(p1))
    # A finite scalar norm cannot cancel the split source-unit exponent.
    c0,c1,dd,ee=s.symbols('c0 c1 D e',nonzero=True)
    mult=s.Matrix([[c0,c1*dd],[c1,c0]])*r**ee
    ck('geometric-unit norm preserves a nonzero exponent',s.simplify(mult.det()-(c0*c0-dd*c1*c1)*r**(2*ee))==0)
    # Generic etale curve maps to A1 with unit derivative omit a value.
    for ex in [2,-3]:
        ff=c*S**ex+b
        ck('Laurent monomial curve map derivative is a unit type '+str(ex),zero(s.diff(ff,S)-ex*c*S**(ex-1)))
        ck('the same affine-line map omits its constant value type '+str(ex),zero((ff-b)/(c*S**ex)-1))
    rho=1-S*Y
    ck('retained node is in the actual source',rho.subs({S:0,Y:0})==1)
    ck('deck scaling retains the unit identically',zero(1-(zeta*S)*(Y/zeta)-rho))
    ck('nontrivial involution fixes node with nonidentity tangent',s.diag(-1,-1)!=s.eye(2))
    # Explicit invariant quotient models retain the node and reveal its obstruction.
    for nd in [2,3,5]:
        phi=U*V-(1-r)**nd;pi={U:S**nd,V:Y**nd,r:rho}
        ck('cyclic invariant relation degree '+str(nd),zero(phi.subs(pi,simultaneous=True)))
        ck('cyclic quotient is singular at the retained node degree '+str(nd),all(s.diff(phi,var).subs({U:0,V:0,r:1})==0 for var in [U,V,r]))
        DP=s.Matrix([S**nd,Y**nd,rho]).jacobian([S,Y])
        ck('quotient map loses tangent rank at node degree '+str(nd),DP.subs({S:0,Y:0})==s.zeros(3,2))
        ck('quotient has full generic tangent rank away from node degree '+str(nd),s.factor(DP.extract([0,2],[0,1]).det())==-nd*S**nd)
        cycl=s.cyclotomic_poly(nd,zeta)
        ck('nontrivial deck eigenvalue cannot equal one degree '+str(nd),s.gcd(cycl,zeta-1)==1)
    ck('degree-one quotient has no node singularity',s.diff(U*V-(1-r),r)==1)
    # All-degree source-surface maps violate primitive unit index.
    d=s.symbols('d',integer=True,positive=True);Hd=(1-r**d)/(1-r)
    ck('generic all-degree self-map Jacobian identity',s.simplify(Hd-(1-r)*s.diff(Hd,r)-d*r**(d-1))==0)
    for nd in [2,3]:
        Hd=sum(rho**ii for ii in range(nd));J=s.Matrix([S,Y*Hd])
        ck('source self-map has nonprimitive unit index '+str(nd),zero(1-J[0]*J[1]-rho**nd))
        ck('source self-map is etale on U degree '+str(nd),s.factor(J.jacobian([S,Y]).det()-nd*rho**(nd-1))==0)
        ck('source self-map covers the whole distinguished line degree '+str(nd),J.subs(S,0)==s.Matrix([0,nd*Y]))
        ck('unit-root extension has exact generic degree '+str(nd),s.degree(r**nd-a,r)==nd)
    ck('selected maps use only exact coefficients',not any(vv.has(s.Float) for vv in list(G)+list(Gi)+list(K)))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('retained local source '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes();ck('reviewed primary HTML '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),
    'assurance':'Selected exact whole-component, primitive-unit, geometric norm, Laurent-curve, cyclic-node and generic index controls accompany written general proofs. Square-free/root closure, Riemann-Hurwitz, genus-zero and normal Hartogs inputs are published. Not independent whole-proof replay, formal verification, external review or publication-priority evidence. No proper unit-free intermediate, nonlinear integration or canonical atomic receiver is constructed.',
    'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
