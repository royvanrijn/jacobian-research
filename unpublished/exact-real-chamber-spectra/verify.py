"""Exact finite controls for the real-chamber paper; not an all-degree proof."""
import json
from pathlib import Path
import sympy as s

W,T=s.symbols("W T")
x,u,gamma,C=s.symbols("x u gamma C")
p,h,pp=s.symbols("p h pp")
assert s.Matrix([[0,gamma,u],[0,0,1],[gamma,0,x]]).det()==gamma**2
assert s.Matrix([[W*pp+gamma,W],[pp,1]]).det()==gamma

def seed(n):
    k=n//2
    roots=([2]+[r for j in range(3,k+2) for r in (j,1-j)] if n%2
           else [3,4]+[r for j in range(5,k+3) for r in (j,1-j)])
    G=s.prod(W-r for r in [0,-1]+roots).expand()
    lam=s.diff(G,W).subs(W,0)-s.diff(G,W).subs(W,1)
    H=(G-s.diff(G,W).subs(W,0)*W)/lam
    return G,lam,s.expand(H)

records=[]
for n in range(3,9):
    G,lam,H=seed(n);p=s.diff(H,W);q=W*p-H;kappa=s.diff(H,W,2).subs(W,1)
    assert lam!=0 and kappa!=-2
    assert [H.subs(W,0),p.subs(W,0),H.subs(W,1),p.subs(W,1)]==[0,0,0,-1]
    a=-(1+kappa)/(2+kappa)
    qbar=s.div(q,W**2,W);assert qbar[1]==0
    assert qbar[0].subs(W,1)==-1
    assert s.diff(qbar[0],W).subs(W,1)==kappa+2
    assert s.simplify(-1+(1+a)*(kappa+2))==0
    b=-s.diff(G,W).subs(W,0)/lam
    for r in s.polys.polytools.ground_roots(G):
        gam=b-p.subs(W,r);xx=1/gam;yy=r-gam;zz=(gam-1-a*xx*yy)/xx**2
        uu=1+xx*yy;gg=1+a*xx*yy+xx**2*zz;ww=uu*gg
        assert s.simplify(ww-r)==0 and s.simplify(gg-gam)==0
        assert s.simplify((1+p.subs(W,ww)/gg)/xx-b)==0
        assert s.simplify((uu+q.subs(W,ww)/gg**2)/xx**2)==0
    counts=None
    if n<=6:
        slope=b+s.Rational(1,10000)
        E=s.Poly(H-slope*W+T,W)
        D=s.Poly(s.discriminant(E.as_expr(),W),T)
        assert s.gcd(D,D.diff()).degree()==0
        intervals=D.intervals(eps=s.Rational(1,100000))
        assert len(intervals)==n-1 and all(mult==1 for _,mult in intervals)
        gaps=[intervals[0][0][0]-1]+[(intervals[j][0][1]+intervals[j+1][0][0])/2 for j in range(len(intervals)-1)]+[intervals[-1][0][1]+1]
        counts=sorted(set(int(s.Poly(E.as_expr().subs(T,t),W).count_roots(-s.oo,s.oo)) for t in gaps))
        assert counts==list(range(n%2,n+1,2)),(n,counts)
    records.append({"degree":n,"rational_all_real_roots":n,"exact_chamber_counts":counts})
report={"scope":"Selected exact seed, reconstruction, divisibility and chamber controls; not the all-degree theorem or a novelty audit.","checks":records}
Path(__file__).with_name("symbolic_checks.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
