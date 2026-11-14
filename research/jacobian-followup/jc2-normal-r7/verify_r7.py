"""Exact two-class nonimmersion classification for the normal F2 r=7 row.

Default includes literal reverse-ideal certificates, in six 30-second workers.
--compact-only checks the small identities and previously pinned receipts,
and explicitly does not claim a new replay of the large local certificate.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import sympy as s

ROOT=Path(__file__).resolve().parent
RESEARCH=ROOT.parents[1]
sys.path.insert(0,str(RESEARCH/'archive/non-elliptic/scripts'))
from verify_f2_degree_6_10_gap import construct_obstructions


def readpoly(text,loc):return s.sympify(text.replace('^','**'),locals=loc)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--compact-only',action='store_true');args=parser.parse_args()
    a,b,c,d,B,D,E,t,z=s.symbols('a b c d B D E t z')
    variables=(a,b,c,d,B,D,E,t,z);loc={str(v):v for v in variables}
    R,rho,_,_,powers=construct_obstructions()
    rho={j:f.as_expr() for j,f in rho.items()}
    assert s.expand((rho[5]+2*a*rho[3]/3+c*rho[1]/3+5*d**3/81).subs(b,0))==0
    p=t**6+a*t**4+b*t**3+c*t**2+d*t
    q=s.expand(sum(co*sum(powers[j][n].as_expr()*t**(2*j-n) for n in range(2*j+1)) for j,co in [(5,1),(4,B),(2,D),(1,E)]))
    assert s.diff(p,t).subs(t,0)==d
    # The tcrit=0 chart: normalize b=1, then eliminate D and E with unit coefficients.
    zero=json.loads((ROOT/'critical-zero-elimination.json').read_text())
    Dzero=readpoly(zero['D'],loc);Ezero=readpoly(zero['E'],loc)
    assert s.expand(s.diff(q,t).subs({t:0,b:1,d:0}).subs(D,Dzero))==0
    assert s.expand(rho[1].subs({b:1,d:0}).subs(D,Dzero).subs(E,Ezero))==0
    fzero=[s.Poly(s.expand(rho[j].subs({b:1,d:0}).subs(D,Dzero).subs(E,Ezero)),B,c,a).clear_denoms()[1].as_expr() for j in [3,5,7,9]]
    assert all(s.expand(f-readpoly(g,loc))==0 for f,g in zip(fzero,zero['equations']))
    multipliers=[readpoly(r.split(' ',2)[2],loc) for r in (ROOT/'critical-zero-reduced.log').read_text().splitlines() if r.startswith('CERT ')]
    assert len(multipliers)==4 and s.expand(sum(f*g for f,g in zip(fzero,multipliers)))==1
    # The tcrit=1 chart: d and E are eliminated with literal unit coefficients.
    data=json.loads((ROOT/'critical_one_elimination.json').read_text())
    dv=readpoly(data['d'],loc);ev=readpoly(data['E'],loc)
    assert s.expand(s.diff(p,t).subs(t,1).subs(d,dv))==0
    assert s.expand(s.diff(q,t).subs(t,1).subs(d,dv).subs(E,ev))==0
    fs=[s.Poly(s.expand(rho[j].subs(d,dv).subs(E,ev)),D,B,c,a,b,z).clear_denoms()[1].as_expr() for j in [1,3,5,7,9]]
    assert all(s.expand(f-readpoly(g,loc))==0 for f,g in zip(fs,data['equations']))
    H=[readpoly(r.split('=',1)[1],loc) for r in (ROOT/'r7-critical-one.log').read_text().splitlines() if r.startswith('H[')]
    assert len(H)==6 and s.expand(H[0]-(2*z+3)**3*(4*z+15))==0
    parameterization={v:s.solve(h,v)[0] for v,h in zip([b,a,c,B,D],H[1:])}
    for f in fs+[b*z-1]:
        assert s.rem(s.expand(f.subs(parameterization,simultaneous=True)),H[0],z)==0
    if not args.compact_only:
        subprocess.run([sys.executable,str(ROOT/'run_lift_checks.py')],check=True,timeout=180)
    meta=json.loads((ROOT/'lift_receipt.json').read_text())
    for column in range(1,7):
        rec=json.loads((ROOT/('lift-checked-%d.json'%column)).read_text())
        assert rec['literal_identity'] and rec['matrix_sha256']==meta['sha256']
        assert rec['generators_sha256']==hashlib.sha256((ROOT/'critical_one_elimination.json').read_bytes()).hexdigest()
        assert rec['targets_sha256']==hashlib.sha256((ROOT/'r7-critical-one.log').read_bytes()).hexdigest()
    cases=[]
    for value in [s.Rational(-3,2),s.Rational(-15,4)]:
        params={v:s.cancel(f.subs(z,value)) for v,f in parameterization.items()}
        params[d]=s.expand(dv.subs(params));params[E]=s.expand(ev.subs(params))
        assert all(s.expand(rho[j].subs(params))==0 for j in [1,3,5,7,9])
        last=s.cancel(rho[11].subs(params));assert last!=0
        pp=s.expand(p.subs(params));qq=s.expand(q.subs(params))
        assert s.degree(pp,t)==6 and s.degree(qq,t)==10
        u=s.Symbol('u')
        field=s.QQ.frac_field(t)
        norm=s.gcd(s.Poly(pp.subs(t,u)-pp,u,domain=field),s.Poly(qq.subs(t,u)-qq,u,domain=field)).monic().as_expr()
        assert norm==u-t
        derivative=s.Poly(s.gcd(s.diff(pp,t),s.diff(qq,t)),t).monic().as_expr()
        assert derivative==t-1
        cross=s.cancel((s.diff(pp,t,2)*s.diff(qq,t,3)-s.diff(qq,t,2)*s.diff(pp,t,3)).subs(t,1))
        assert s.diff(pp,t,2).subs(t,1)!=0 and cross!=0
        fiber=s.Poly(s.gcd(pp-pp.subs(t,1),qq-qq.subs(t,1)),t).monic().as_expr()
        assert s.expand(fiber-(t-1)**2)==0
        dp=s.div(pp.subs(t,u)-pp,u-t,u)[0];dq=s.div(qq.subs(t,u)-qq,u-t,u)[0]
        conductor=s.Poly(s.resultant(dp,dq,u),t).clear_denoms()[1].primitive()[1].as_expr()
        assert s.degree(conductor,t)==24
        if value==s.Rational(-3,2):
            multiple=3*t*t+3*t+1
            remainder,rem=s.div(conductor,(t-1)**2*multiple**3,t);assert rem==0
            assert s.degree(remainder,t)==16 and s.degree(s.gcd(remainder,s.diff(remainder,t)),t)==0
            assert s.degree(s.gcd(remainder,(t-1)*multiple),t)==0
            assert s.rem(pp,multiple,t)==-s.Rational(1,9)
            assert s.rem(qq,multiple,t)==s.Rational(49,59049)
            otherfiber=s.Poly(s.gcd(pp+s.Rational(1,9),qq-s.Rational(49,59049)),t).monic().as_expr()
            assert s.expand(otherfiber-multiple/3)==0
            singularities='A2 + A5 + 8 A1'
        else:
            remainder,rem=s.div(conductor,(t-1)**2,t);assert rem==0
            assert s.degree(remainder,t)==22 and s.degree(s.gcd(remainder,s.diff(remainder,t)),t)==0
            assert s.degree(s.gcd(remainder,t-1),t)==0
            singularities='A2 + 11 A1'
        cases.append({'inverse_b':str(value),'parameters':{str(v):str(f) for v,f in params.items()},
                      'p':str(pp),'q':str(qq),'rho11':str(last),'normalization_gcd':str(norm),
                      'derivative_gcd':str(derivative),'cusp_jet_determinant':str(cross),
                      'conductor':str(conductor),'affine_singularities':singularities})
    # Three distinct smooth points with one vertical tangent would make
    # p-P0 the square of a monic cubic. Centering and b=1 force this family.
    U=s.Symbol('u');vertical=json.loads((ROOT/'vertical-triple-input.json').read_text())
    vsub={a:2*U,b:1,c:U**2,d:U}
    cubic=t**3+U*t+s.Rational(1,2)
    assert s.expand(p.subs(vsub)-(cubic**2-s.Rational(1,4)))==0
    qv=s.expand(q.subs(vsub));qrem=s.Poly(s.rem(qv,cubic,t),t)
    expressions=[s.expand(rho[j].subs(vsub)) for j in [1,3,5,7,9]]+[qrem.coeff_monomial(t),qrem.coeff_monomial(t*t)]
    vf=[s.Poly(f,E,D,B,U).clear_denoms()[1].as_expr() for f in expressions]
    vloc={**loc,'u':U}
    assert all(s.expand(f-readpoly(g,vloc))==0 for f,g in zip(vf,vertical['equations']))
    vm=[readpoly(r.split(' ',2)[2],vloc) for r in (ROOT/'vertical-triple-lift.log').read_text().splitlines() if r.startswith('CERT ')]
    assert len(vm)==7 and s.expand(sum(f*g for f,g in zip(vf,vm)))==1
    subprocess.run([sys.executable,str(ROOT/'local_excess.py')],check=True,timeout=10)
    local=json.loads((ROOT/'local-excess.json').read_text())
    assert local['required_Orevkov_excess']==3 and local['A2_plus_11_nodes_excluded_under_complete_one_row_ledger']
    pi,w=s.symbols('pi w');aa=s.Symbol('aa');zz=w**3+aa*w
    assert s.expand(-4*aa**3-27*zz**2+(3*w**2+aa)**2*(3*w**2+4*aa))==0
    result={'scope':'Birational polynomial (6,10) pairs with first odd infinity gap 21, in the archived normalization; classification of affine nonimmersion only',
            'scheme_equations':[str(f) for f in H],'scheme_length':4,
            'support_equation':'(2*z+3)*(4*z+15)=0; z=1/b after centering the marked critical parameter at 1 by scaling',
            'critical_zero_unit_identity':True,'both_ideal_inclusions_verified':True,
            'vertical_triple_tangent_unit_identity':True,
            'case_II_excluded_conditionally_on_complete_one_row_ledger':True,
            'case_I_required_excess_split_under_same_assumption':[1,2],
            'ramified_and_unramified_companion_factorization':True,
            'large_reverse_certificate_replayed_this_run':not args.compact_only,'large_certificate_sha256':meta['sha256'],
            'cases':cases,'boundary':'No full parameter-family topology, degree-six exclusion, (75,125) exclusion, coordinate-degree improvement or JC2 proof.'}
    (ROOT/'classification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS EXACT TWO-CLASS R7 NONIMMERSION CLASSIFICATION')
    print('Case singularities:',[row['affine_singularities'] for row in cases])


if __name__=='__main__':main()
