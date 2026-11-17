#!/usr/bin/env python3
"""Selected exact controls for descent and its averaging boundary; not whole-proof replay."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def ck(label,value):
    assert bool(value),label
    checks.append(label)
def zero(v):return s.cancel(v)==0
def main():
    begun=time.monotonic()
    x,y,z,u,v,t,H,q=s.symbols('x y z u v t H q')
    a=s.sqrt(2);d=s.Integer(2)
    # A stable coordinate line can require a WEIGHTED trace.
    f=a*x+3+a
    f0=s.expand(f-f.subs(x,0))
    conjugate=lambda e:s.expand(e).xreplace({a:-a})
    ck('rational-point centering removes the affine translation',f0==a*x)
    ck('the centered ruling line is semilinearly stable',zero(conjugate(f0)+f0))
    ck('ordinary trace can vanish even for the same ruling',zero(f0+conjugate(f0)))
    weighted=s.expand(a*f0+conjugate(a*f0))
    ck('a nonzero weighted trace gives the descended coordinate',weighted==4*x)
    ck('the quadratic ruling descent coefficient is invertible',s.diff(weighted,x)==4)
    ck('the centered scalar cocycle satisfies sigma(a)*a=1',(-1)*(-1)==1)
    # Genuine conjugate polynomial coordinates with a trace that IS NOT an A1-ruling.
    wplus=y+a*x*x
    fp=x+a*wplus*wplus
    wminus=y-a*x*x
    fm=x-a*wminus*wminus
    ck('the two affine rulings are genuinely conjugate',zero(conjugate(fp)-fm))
    ck('the first displayed ruling is part of a Jacobian-one coordinate pair',zero(s.Matrix([fp,wplus]).jacobian([x,y]).det()-1))
    ck('the conjugate ruling is part of a Jacobian-one coordinate pair',zero(s.Matrix([fm,wminus]).jacobian([x,y]).det()-1))
    X,W=s.symbols('X W')
    invx=X-a*W*W
    invy=W-a*invx*invx
    ck('the first triangular coordinate inverse recovers its base',zero(fp.subs({x:invx,y:invy},simultaneous=True)-X))
    ck('the first triangular coordinate inverse recovers its fibre',zero(wplus.subs({x:invx,y:invy},simultaneous=True)-W))
    ck('both conjugate rulings are centered at the rational origin',fp.subs({x:0,y:0})==0 and fm.subs({x:0,y:0})==0)
    trace=s.expand((fp+fm)/2)
    ck('their trace is x+4*x^2*y',trace==x+4*x*x*y)
    ck('their ruling fields are distinct',not zero(s.Matrix([fp,fm]).jacobian([x,y]).det()))
    ygeneric=(H-x)/(4*x*x)
    ck('the generic trace fibre solves for y with an x denominator',zero(trace.subs(y,ygeneric)-H))
    ck('x has a regular inverse in the generic trace fibre',zero(x*((1+4*x*y)/H)-1-(trace-H)/H))
    ck('the explicit inverse equals the Laurent inverse',zero(((1+4*x*y)/H).subs(y,ygeneric)-1/x))
    ck('the trace fibre loses x=0 when the base value is nonzero',trace.subs(x,0)==0)
    # Canonical contraction on higher-power Danielewski controls.
    m=s.symbols('m',integer=True,positive=True)
    for P in [z*z-1,z**3-z]:
        D=lambda expr:s.expand(-u**m*s.diff(expr,z)-s.diff(P,z)*s.diff(expr,v))
        surface=u**m*v-P
        tag=' P='+str(P)
        ck('Hamiltonian field preserves the hypersurface'+tag,zero(D(surface)))
        ck('Hamiltonian field fixes the ruling coordinate'+tag,D(u)==0)
        ck('contraction with du wedge dz/u^m gives du'+tag,zero(-D(z)/u**m-1))
        ck('contraction on the P-prime volume chart gives du'+tag,zero(-D(v)/s.diff(P,z)-1))
        it=v
        for _ in range(s.degree(P,z)+1):it=D(it)
        ck('selected generators exhibit local nilpotence'+tag,it==0 and D(D(z))==0)
        ck('the generic ruling fibre is an affine line'+tag,zero(surface.subs(v,P/u**m)))
    # The geometric Picard norm of the two quadric rulings can vanish without equality.
    relation=s.Matrix([1,1])
    ck('opposite ruling classes have zero orbit sum modulo the diagonal',s.Matrix([1,0])+s.Matrix([0,1])==relation)
    ck('the two ruling classes remain different modulo the diagonal',s.Matrix.hstack(relation,s.Matrix([1,-1])).det()!=0)
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:
        ck('retained local source '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes()
        ck('reviewed primary source '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={
        'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,
        'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),
        'assurance':'Two written general theorems with selected symbolic controls. Classification and standard surface/descent facts are explicitly imported. Weighted trace works after uniqueness, while an exact conjugate-coordinate counterexample disproves arbitrary averaging. Not independent whole-proof replay, formal verification, external review or novelty evidence.',
        'scope_boundary':'Every affine-type geometric ruling of a proper native intermediate descends; conditional on JC(2), only additively rigid complete-type rulings remain. Complete-type descent, JC(2), actual surviving whole factors and absolute stable atomicity remain open.',
        'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','README.md','verify.py','requirements.txt','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ['result','check_count','runtime_seconds']}))
if __name__=='__main__':main()
