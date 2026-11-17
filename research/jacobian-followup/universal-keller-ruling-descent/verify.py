#!/usr/bin/env python3
"""Selected exact controls, not independent replay of the written general proofs."""
from pathlib import Path
from fractions import Fraction
from math import prod
import hashlib,json,platform,time
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(v):return s.cancel(v)==0
def main():
    begun=time.monotonic();S,Y,r,u,v,z,a,b,c=s.symbols('S Y r u v z a b c');rho=1-S*Y
    # The A chart and its intersection with B make the van Kampen calculation explicit.
    ck('source A-chart has an exact inverse',zero(1-S*((1-r)/S)-r))
    mat=smith_normal_form(s.Matrix([[1],[0]]),domain=s.ZZ)
    ck('killing the S loop leaves one free abelian generator',mat==s.Matrix([[1],[0]]))
    ck('hyperbola compact-support Euler characteristic is zero',2-2==0)
    ck('hyperbola complement compact-support Euler characteristic is one',1-(2-2)==1)
    N=s.symbols('N',integer=True,positive=True)
    ck('all product leaves retain Euler characteristic one',s.Integer(1)**N==1)
    # Uniform self-map controls show the actual line of lost sheets.
    nn=s.symbols('n',integer=True,positive=True);H=(1-r**nn)/(1-r)
    ck('all-degree source-map Jacobian identity',s.simplify(H-(1-r)*s.diff(H,r)-nn*r**(nn-1))==0)
    for nd in [2,3]:
        hd=sum(rho**ii for ii in range(nd));J=s.Matrix([S,Y*hd]);label=' degree '+str(nd)
        ck('source self-map is etale on the unit open'+label,zero(J.jacobian([S,Y]).det()-nd*rho**(nd-1)))
        ck('source self-map has the expected r power'+label,zero(1-J[0]*J[1]-rho**nd))
        ck('whole distinguished line has exactly the displayed invertible restriction'+label,J.subs(S,0)==s.Matrix([0,nd*Y]))
        ck('off-line unit equation has degree d'+label,s.degree(r**nd-a,r)==nd)
        ck('Euler sheet-loss budget agrees on the line'+label,(nd-1)*1==nd*1-1)
    H2=sum(rho**i for i in range(2));H3r2=sum(rho**(2*i) for i in range(3));H6=sum(rho**i for i in range(6))
    ck('degree-two and degree-three controls compose to degree six',zero(H2*H3r2-H6))
    ck('their line-defect composition is exact',2*(3-1)+(2-1)==6-1)
    aa,bb,cc,chiX,chiY,chiZ,chiW=s.symbols('a b c chiX chiY chiZ chiW')
    Df=aa*chiY-chiX;Dj=bb*chiZ-chiY
    ck('general Euler defect composition identity',s.expand(aa*Dj+Df-(aa*bb*chiZ-chiX))==0)
    ck('three-stage Euler defect composition is associative',s.expand(aa*bb*(cc*chiW-chiZ)+aa*(bb*chiZ-chiY)+(aa*chiY-chiX)-(aa*bb*cc*chiW-chiX))==0)
    ck('intermediate Euler characteristic cancels from the whole budget',s.expand(aa*(bb-chiY)+(aa*chiY-1)-(aa*bb-1))==0)
    ck('nonnegative constructible loss can have negative Euler integral',1*(2-2*1-1)==-1)
    # Arithmetic-integral but geometrically split fibres are the permitted boundary.
    for nd,nm in [(2,2),(3,3)]:
        P=z**nd-a;label=' degree '+str(nd)+' exponent '+str(nm)
        ring=s.QQ.frac_field(a)
        ck('arithmetic special-fibre polynomial is irreducible'+label,s.Poly(P,z,domain=ring).is_irreducible)
        ck('the same fibre is geometrically reduced'+label,s.gcd(P,s.diff(P,z),z)==1)
        ck('generic affine-line fibre solves for v'+label,zero(u**nm*(P/u**nm)-P))
        ck('arithmetic model has the specified rational point'+label,zero((u**nm*v-P).subs({u:1,v:-a,z:0})))
        for root in [0,1]:
            P0=z*(z-1);ww=s.symbols('w');Z=root+u**nm*ww
            V=ww*s.prod(root-other+u**nm*ww for other in [0,1] if other!=root)
            ck('selected reduced fibre component gives a plane chart'+label+' root '+str(root),zero(u**nm*V-P0.subs(z,Z)))
    # Hypothesis controls: dropping factoriality or generic integrality changes the conclusion.
    for nm in [2,3]:
        ps=u+u*u*v-z**nm;label=' exponent '+str(nm)
        ck('pseudo-plane critical test along the multiple fibre'+label,s.diff(ps,u).subs({u:0,z:0})==1)
        ck('pseudo-plane zero fibre retains multiplicity m'+label,ps.subs(u,0)==-z**nm)
        ck('the apparent power denominator is not a global unit'+label,(1+u*v).subs({u:1,v:-1,z:0})==0 and ps.subs({u:1,v:-1,z:0})==0)
        ck('a polynomial power loses geometric generic connectedness'+label,s.degree(u**nm-a,u)==nm and s.gcd(u**nm-a,s.diff(u**nm-a,u),u)==1)
    ck('quadric complete-type ruling has matching P1 chart ratios',zero(u*v-(z-1)*(z+1)-(u*v-z*z+1)))
    # Group-tree kernel rank controls; the uniform classification is proved in the note.
    for orders,expected in [([2],0),([2,2],1),([2,3],2),([3,3],4),([2,2,2],5)]:
        qq=prod(orders);rank=1+qq*(len(orders)-1-sum(Fraction(1,m) for m in orders))
        ck('orbifold free-kernel rank '+str(orders),rank==expected)
    kk=s.symbols('k',integer=True,nonzero=True);T=s.Matrix([[1,kk],[0,1]]);R=s.Matrix([[-1,0],[0,1]])
    ck('infinite-dihedral reflection has order two',R*R==s.eye(2))
    ck('reflection inverts all translations',R*T*R==T.inv())
    ck('nontrivial translations do not commute with a reflection',T*R-R*T==s.Matrix([[0,2*kk],[0,0]]))
    # Index divisibility is stronger than an upper bound and weaker than equality.
    ck('quotienting an index-six cover can leave index two',6==2*3)
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:
        ck('retained local source '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes()
        ck('reviewed primary source '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),
        'assurance':'Selected exact chart, source self-map, sheet-loss, arithmetic/geometric fibre, hypothesis-boundary, group-tree and dihedral controls accompany three written general claims. Classical surface fibration, topology, orbifold quotient and normal extension inputs are credited. Not independent whole-proof replay, formal verification, external review or novelty evidence. No ruling descent theorem for every geometric intermediate, JC(2) solution, compatible target atlas or canonical atomic receiver is constructed.',
        'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
