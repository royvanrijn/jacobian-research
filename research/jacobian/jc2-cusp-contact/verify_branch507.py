#!/usr/bin/env python3
"""Exact finite branch507 certificates; no historic F2 DAG reconstruction."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import platform
import time
import sympy as S

H=Path(__file__).resolve().parent
R=H.parents[1]
checks=[]
def check(name,c):
    if not c: raise AssertionError(name)
    checks.append(name)
def zero(name,expr):
    check(name,S.cancel(S.expand(expr))==0)

def hull(points):
    ps=sorted(set(points))
    def cross(o,a,b):
        return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lo=[]
    for p in ps:
        while len(lo)>1 and cross(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    up=[]
    for p in reversed(ps):
        while len(up)>1 and cross(up[-2],up[-1],p)<=0: up.pop()
        up.append(p)
    return lo[:-1]+up[:-1]

def area(p):
    return Fraction(abs(sum(x*v-y*u for (x,y),(u,v)
                            in zip(p,p[1:]+p[:1]))),2)
def mixed(a,b):
    c=hull([(x+u,y+v) for x,y in a for u,v in b])
    return area(c)-area(a)-area(b)

def main():
    start=time.monotonic()
    manifest=json.loads((H/'sources/bkk-manifest.json').read_text())
    check('pinned primary BKK source',
          sha256((H/'sources'/manifest['file']).read_bytes()).hexdigest()==manifest['sha256'])
    P=hull([(i,j) for i in range(76) for j in range(76-i) if 5*i-j<=15])
    Q=hull([(i,j) for i in range(126) for j in range(126-i) if 5*i-j<=25])
    check('exact P envelope',P==[(0,0),(3,0),(15,60),(0,75)])
    check('exact Q envelope',Q==[(0,0),(5,0),(25,100),(0,125)])
    check('direct generic degree upper bound2175',mixed(P,Q)==2175)
    expected=[1347,1384,1554,1729,1899,1999,2009,2019]
    rows=[]
    for L,upper in zip((39,40,45,50,55,60,65,70),expected):
        T=min(375,13*L//2)
        W=hull([(i,j) for i in range(376) for j in range(376-i)
                if 5*i-j<=L and 5*i+2*j<=13*L])
        mv=mixed(P,W)
        check('integer-hull mixed volume L='+str(L),mv==12*L+15*T)
        check('integer-hull facet supports L='+str(L),
              max(5*i-j for i,j in W)==L and max(i+j for i,j in W)==T)
        check('terminal local intersection budget L='+str(L),
              Fraction(72*(5*L-15),60)==6*L-18)
        check('refined sheet-degree bound L='+str(L),
              (mv-(6*L-12))/3==upper==2*L+5*T+4)
        rows.append(dict(L=L,contact=13*L,total_cap=T,
                         W_hull=[list(p) for p in W],mixed_volume=int(mv),
                         eta_intersection_lower_bound=6*L-12,
                         raw_degree_bound=int(mv/3),refined_degree_bound=upper))
    check('uniform F2 geometric-degree bound2019',max(expected)==2019)
    # Endpoint height transport, as actual monomial equations in the retained chart.
    for aa,bb in ((4,1),(21,6)):
        ell=aa-bb
        check('P terminal endpoint local weight '+str((aa,bb)),
              5*(15-ell)+12*bb==72)
    for aa,bb in ((1,0),(18,5),(35,10)):
        ell=aa-bb
        check('g terminal endpoint local weight '+str((aa,bb)),
              5*(25-ell)+12*bb==120)
    ss=S.symbols('ss')
    terminal=S.expand(ss*(1+ss)**5-(ss**2+S.Rational(5,3)*ss+S.Rational(5,9))**3)
    check('actual terminal W polynomial has degree3',S.degree(terminal,ss)==3)
    zero('terminal W has no common fast-root with P',terminal.subs(ss,-1)-S.Rational(1,729))
    zero('terminal W nonzero slow-root coefficient',terminal.subs(ss,0)+S.Rational(125,729))
    for L in (39,40,45,50,55,60,65,70):
        check('terminal W cleared monomials have stated local weight '+str(L),
              all(L-3-12*j>=0 and 5*(L-3-12*j)+12*(5*j)==5*L-15
                  for (j,),coef in S.Poly(terminal,ss).terms()))
    q,z=S.symbols('q z')
    resultant=S.resultant(z**6+q**15,z**15+q**36,z)
    zero('local weighted bound control resultant',resultant-q**216*(1+q**3)**3)

    t,a,b,c,d,mu,eta,nu=S.symbols('t a b c d mu eta nu')
    pp=t**3+a*t+mu
    qq=t**5+b*t**4+eta*t**3+c*t**2+(d+eta*a)*t+nu
    difference=S.Poly(pp**5-qq**3,t)
    coeff=[difference.coeff_monomial(t**j) for j in range(14,7,-1)]
    G=S.groebner(coeff,b,eta,c,d,nu,mu,a)
    basis=[b,3*eta-5*a,3*c-5*mu,9*d+10*a*a,
           9*nu-10*a*mu,9*mu*mu-a**3,a*a*mu,a**5]
    check('zero seven-jet full coefficient ideal',
          list(G)==basis)
    triangular={b:0,eta:5*a/3,c:5*mu/3,d:-10*a*a/9,nu:10*a*mu/9}
    for i in range(5):
        zero('seven-jet triangular pivot '+str(i),coeff[i].subs(triangular))
    zero('first residual raw jet coefficient',
         coeff[5].subs(triangular)+S.Rational(5,27)*(a**3-9*mu**2))
    zero('second residual raw jet coefficient',
         coeff[6].subs(triangular)+S.Rational(5,9)*a*a*mu)
    zero('raw jet nilpotence identity',
         a**5-9*(mu*a*a*mu-a*a*(mu*mu-a**3/9)))
    standard=list(map(S.sympify,(1,a,a*a,a**3,a**4,mu,a*mu)))
    check('raw jet length7 standard monomials',
          all(G.reduce(x)[1]==x for x in standard))
    # Directly eliminate using an order placing transports first, then compare.
    all_elim=S.groebner(coeff,eta,nu,mu,b,d,c,a)
    target_generators=[p for p in all_elim if not p.has(mu,eta,nu)]
    target=S.groebner(target_generators,b,d,c,a,domain=S.QQ)
    expected_target=[b,d+S.Rational(10,9)*a*a,
                     c*c-S.Rational(25,81)*a**3,a*a*c,a**5]
    check('zero-jet elimination equals exact Euler eligibility scheme',
          list(target)==expected_target)
    reverse={mu:3*c/5,eta:5*a/3,nu:2*a*c/3}
    for i,poly in enumerate(coeff):
        check('explicit reverse scheme map coefficient '+str(i),
              target.reduce(S.expand(poly.subs(reverse)))[1]==0)
    check('no geometric non-E8 zero-seven-jet row',
          G.reduce(a**5)[1]==0 and G.reduce(mu**3)[1]==0)

    receipt=dict(as_of='2026-11-14',outcome='PASS',checks=checks,check_count=len(checks),
                 geometric_degree_atlas=rows,
                 zero_seven_jet_basis=list(map(str,G)),
                 elimination_Euler_basis=list(map(str,target)),
                 zero_jet_scheme_length=7,
                 files={n:sha256((H/n).read_bytes()).hexdigest()
                        for n in ('BRANCH507.md','verify_branch507.py',
                                  'sources/bkk-manifest.json')},
                 source_hashes={n:sha256((R/n).read_bytes()).hexdigest() for n in (
                    'foundations/plane-jc/F2_AFFINE_K1_KELLER_PULLBACK.md',
                    'foundations/plane-jc/F2_AFFINE_K1_CARRIER_JET_FACTORIZATION.md',
                    'foundations/plane-jc/F2_AFFINE_PURITY_PUNCTURE_ATTACHMENT.md',
                    'foundations/plane-jc/F2_75_125_DERIVATION.md')},
                 mathematical_boundary='BKK and toric/local intersection arguments are written imports/proofs. Zero-jet matching does not establish global attachment. No full F2 closure.',
                 python=platform.python_version(),sympy=S.__version__,
                 elapsed_seconds=round(time.monotonic()-start,4))
    (H/'branch507-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',checks=len(checks),
                         geometric_degree_bound=2019,contact507_bound=1347,
                         elapsed_seconds=receipt['elapsed_seconds'])))

if __name__=='__main__': main()
