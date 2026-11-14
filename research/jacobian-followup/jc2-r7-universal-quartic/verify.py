"""Universal quartic contact identities and two sparse Fraction unit replays.

No Gröbner basis computation is part of this checker. The grading proof
for exact coefficient filtering and inherited geometry are written input.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import time

import sympy as S
from sympy.polys.rings import ring
from weighted_contact import (P,a,b,c,d,B,D,E,u,v,w,coefficient,rho,
                              columns,degrees,weights,remainders,pp,qq)
from contact import contact_matrix,selected_minor,z_coefficient,R as FULL,u as FULL_U

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(name,value):
    if not value:raise AssertionError(name)
    CHECKS.append(name)


def zero(name,value):check(name,S.expand(value)==0)


def decode(rows,n):
    result={}
    for exponent,coefficient in rows:
        m=tuple(exponent);q=Fraction(coefficient)
        assert len(m)==n and all(type(e) is int and e>=0 for e in m)
        assert m not in result and q
        result[m]=q
    return result


def multiply(f,g):
    result={}
    for m,a0 in f.items():
        for n,b0 in g.items():
            key=tuple(x+y for x,y in zip(m,n))
            result[key]=result.get(key,Fraction(0))+a0*b0
    return {m:q for m,q in result.items() if q}


def add(f,g):
    result=dict(f)
    for m,q in g.items():result[m]=result.get(m,Fraction(0))+q
    return {m:q for m,q in result.items() if q}


def fraction_unit(generators,multipliers,n):
    total={}
    for f,g in zip(generators,multipliers):
        total=add(total,multiply(decode(f,n),decode(g,n)))
    return total


def primitive(expression,variables):
    return S.Poly(S.expand(expression),*variables).clear_denoms()[1].primitive()[1].as_expr()


def main():
    start=time.monotonic()
    sa,sb,sc,sd,sB,sD,sE,su,sv,sw=S.symbols('a b c d B D E u v w')
    check('column grading reconstructed',all(
        all(sum(x*y for x,y in zip(weights,m))==deg-n for m in coefficient_poly)
        for deg,column in zip(degrees,columns) for n,coefficient_poly in enumerate(column)))
    check('doubled-quartic remainder grading',all(
        sum(x*y for x,y in zip((1,2,3,4),m))+r==n
        for n,rows in enumerate(remainders) for r,poly in enumerate(rows) for m in poly))
    zero('gap parity excludes the b0,d-nonzero chart',
         (rho[5].as_expr()+2*sa*rho[3].as_expr()/3+sc*rho[1].as_expr()/3).subs(sb,0)+5*sd**3/81)
    zero('u cubic coefficient universally forced zero',
         coefficient((0,2,3,4),10).as_expr()+80*su**5)
    base={u:P.zero,b:P.one}
    L=-21*sB-20*sa+30*sv;K=-21*sB+20*sa-60*sd
    F10=coefficient((0,1,2,3),10,base).as_expr()
    zero('complete leading parameter split',F10+4*(3*sw-2)**2*L/9)
    wbranch={**base,w:P.ground_new(S.Rational(2,3))}
    F9=coefficient((0,1,2,3),9,wbranch).as_expr()
    zero('L-nonzero branch forces v',F9-4*L*(-4*sa+3*sv+6*sd)**2/81)
    wbranch[v]=4*a/3-2*d
    F8=coefficient((0,1,2,3),8,wbranch).as_expr()
    zero('L-nonzero branch forces c',F8+16*K*(2*sa*sd-sc-3*sd**2)**2/81)
    wbranch[c]=2*a*d-3*d*d
    F7=coefficient((0,1,2,3),7,wbranch).as_expr()
    r1=S.expand(rho[1].as_expr().subs(sb,1))
    zero('L-nonzero branch contradicts rho1',
         F7-S.Rational(8,27)*r1.subs(sc,2*sa*sd-3*sd**2)-4*K/6561)

    Lbranch={**base,v:2*a/3+7*B/10}
    Rexpr=(-3321*sB**3-6300*sB*sa**2+18900*sB*sc+32400*sD
           +4000*sa**3-18000*sa*sc-6000-22500*sw+27000*sw**2)
    F9=coefficient((0,1,2,3),9,Lbranch).as_expr()
    zero('remaining R split',F9-(3*sw-2)**2*Rexpr/6075)
    zero('rho1 E coefficient is a unit',S.diff(r1,sE)-S.Rational(1,3))
    Ebase=-3*r1.subs(sE,0)
    zero('E elimination exact',r1.subs(sE,Ebase))
    zero('R D coefficient is a unit',S.diff(Rexpr,sD)-32400)
    DR=-Rexpr.subs(sD,0)/32400
    ER=Ebase.subs(sD,DR)
    zero('R elimination exact',Rexpr.subs(sD,DR))
    remaining={**Lbranch,D:P.from_expr(DR),E:P.from_expr(ER)}
    other=coefficient((0,1,3,5),8,remaining).as_expr()
    zero('additional minor forces w2/3',other+S.Rational(40,81)*(3*sw-2)**4)
    intersection={**Lbranch,w:P.ground_new(S.Rational(2,3))}
    F8=coefficient((0,1,2,3),8,intersection).as_expr()
    R0=Rexpr.subs(sw,S.Rational(2,3))
    zero('two complete final branches',F8+4*K*r1/15+K**2*R0/5467500)

    certs=json.loads((HERE/'certificates.json').read_text())['units']
    check('both final branches frozen exactly once',len(certs)==2 and {x['branch'] for x in certs}=={'K_zero','K_nonzero'})
    reports=[]
    for cert in certs:
        variables=S.symbols(' '.join(cert['variables']))
        RR,*_=ring(','.join(cert['variables']),S.QQ)
        if cert['branch']=='K_zero':
            BV=(20*sa-60*sd)/21;EV=Ebase.subs(sB,BV)
            sub={**base,w:P.ground_new(S.Rational(2,3)),v:4*a/3-2*d,
                 B:P.from_expr(BV),E:P.from_expr(EV)}
            gaps=[rho[j].as_expr().subs(sb,1).subs(sB,BV).subs(sE,EV) for j in (3,5,7,9)]
            contacts=[coefficient((0,1,2,3),j,sub).as_expr() for j in (7,6)]
        else:
            DV=DR.subs(sw,S.Rational(2,3));EV=Ebase.subs(sD,DV)
            sub={**intersection,D:P.from_expr(DV),E:P.from_expr(EV)}
            gaps=[rho[j].as_expr().subs(sb,1).subs(sD,DV).subs(sE,EV) for j in (3,5,7,9)]
            contacts=[coefficient((0,1,2,3),7,sub).as_expr()]
        generators=[primitive(f,variables) for f in gaps+contacts]
        check(cert['branch']+' generator count',len(generators)==len(cert['generators'])==len(cert['multipliers']))
        for i,(f,encoded) in enumerate(zip(generators,cert['generators'])):
            actual=RR.from_expr(f)
            expected=RR.from_dict({tuple(m):S.QQ(c) for m,c in encoded})
            check(cert['branch']+f' exact generator {i}',actual==expected)
        value=fraction_unit(cert['generators'],cert['multipliers'],len(variables))
        check(cert['branch']+' independent sparse Fraction unit',value=={(0,)*len(variables):Fraction(1)})
        mutated=json.loads(json.dumps(cert['multipliers']))
        first=decode(mutated[0],len(variables));key=(0,)*len(variables)
        first[key]=first.get(key,Fraction(0))+1
        mutated[0]=[[list(m),str(q)] for m,q in sorted(first.items()) if q]
        check(cert['branch']+' mutated unit rejected',
              fraction_unit(cert['generators'],mutated,len(variables))!={(0,)*len(variables):Fraction(1)})
        reports.append(dict(branch=cert['branch'],variables=cert['variables'],generator_count=len(generators),
                            multiplier_terms=[len(f) for f in cert['multipliers']],unit_identity=True))

    # Exact symbolic implementation controls against the earlier full matrix.
    classes=json.loads((HERE.parent/'jc2-normal-r7/classification.json').read_text())['cases']
    generic=coefficient((0,1,2,3),9,{u:P.zero}).as_expr()
    for case,row in enumerate(classes,1):
        numeric={S.Symbol(k):S.Rational(value) for k,value in row['parameters'].items()}
        matrix,_=contact_matrix(S.sympify(row['p']),S.sympify(row['q']))
        atzero=[[f.compose(FULL_U,FULL.zero) for f in row] for row in matrix]
        actual=z_coefficient(selected_minor(atzero,[0,1,2,3]),9).as_expr()
        zero(f'class{case} full-matrix specialization control',generic.subs(numeric)-actual)
    inputs=[HERE/'PROOF.md',HERE/'README.md',HERE/'REVIEW.md',HERE/'certificates.json',HERE/'weighted_contact.py',Path(__file__).resolve(),
            ROOT/'archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py',
            HERE.parent/'jc2-r7-quartic-contact/contact.py',
            HERE.parent/'jc2-r7-quartic-contact/PROOF.md',
            HERE.parent/'jc2-normal-r7/classification.json']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    receipt=dict(schema='universal-r7-quartic-v1',passed=True,check_count=len(CHECKS),checks=CHECKS,
                 elapsed_seconds=round(time.monotonic()-start,3),input_sha256=hashes,units=reports,
                 universal_parameter_forcing=dict(b=1,u=0,v='2*a/3+7*B/10',w='2/3'),
                 universal_quartic_contact='EXCLUDED for every birational degree(6,10) gap21 pair',
                 remaining_single_ramified_row_case='reciprocal cubic with a finite multiple root of p_prime',
                 assurance='Exact universal identities and separate Fraction units. Weighted-filter justification and inherited geometry are written proofs; no formalization or external review. Capped cofactor/basis attempts are not proof inputs.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} UNIVERSAL QUARTIC CHECKS ({receipt["elapsed_seconds"]}s)')
    print('All gap21 quartic conic-contact parameters excluded; two multivariate Fraction units replayed.')


if __name__=='__main__':main()
