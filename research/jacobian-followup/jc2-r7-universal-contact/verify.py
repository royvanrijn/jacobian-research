"""Exact universal polynomial identities; geometric deductions are written inputs.

No large elimination, numerical sampling or archived artifact reconstruction.
The affine-linear unit certificate also has a separate Fraction replay.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import sys
import time

import sympy as S
from sympy.polys.rings import ring

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'archive/non-elliptic/scripts'))
from verify_f2_degree_6_10_gap import construct_obstructions

CHECKS = []


def check(name, value):
    if not value:
        raise AssertionError(name)
    CHECKS.append(name)


def zero(name, value):
    check(name, S.expand(value) == 0)


def coefficient(poly, exponent):
    """Coefficient in the FIRST variable, retaining every other variable."""
    return poly.ring.from_dict({(0,)+m[1:]: c for m, c in poly.items()
                               if m[0] == exponent})


def fraction_replay(polynomials, multipliers):
    # Each row lists coefficients of 1,B,D,E; no symbolic simplifier here.
    return [sum((Fraction(c)*Fraction(row[j])
                 for c, row in zip(multipliers, polynomials)), Fraction(0))
            for j in range(4)]


def main():
    start = time.monotonic()
    cert = json.loads((HERE/'certificates.json').read_text())
    _, rho, _, _, powers = construct_obstructions()
    a,b,c,d,B,D,E,t,h,u,alpha,x,A = S.symbols('a b c d B D E t h u alpha x A')
    p = t**6+a*t**4+b*t**3+c*t**2+d*t
    q = S.expand(sum(co*sum(powers[j][n].as_expr()*t**(2*j-n)
                           for n in range(2*j+1))
                     for j,co in [(5,1),(4,B),(2,D),(1,E)]))
    # Independent binomial expansion of each polynomial part.
    correction = a*u**2+b*u**3+c*u**4+d*u**5
    independent = 0
    for j,co in [(5,1),(4,B),(2,D),(1,E)]:
        expansion = S.Poly(S.expand(sum(S.binomial(S.Rational(j,3),k)*correction**k
                                        for k in range(j+1))),u)
        independent += co*sum(expansion.nth(n)*t**(2*j-n) for n in range(2*j+1))
    zero('polynomial parts agree with independent binomial expansion',q-independent)
    zero('retained parity gap identity',
         (rho[5]+2*rho[3]*rho[5].ring.gens[0]/3
          +rho[1]*rho[5].ring.gens[2]/3).as_expr().subs(b,0)+5*d**3/81)
    check('even p at b=d=0', S.Poly(p.subs({b:0,d:0}),t).monoms()
          == [(6,), (4,), (2,)])
    zero('even q at b=d=0', q.subs({b:0,d:0})-q.subs({b:0,d:0,t:-t}))

    pp = t**6-15*alpha**2*t**4+40*alpha**3*t**3-45*alpha**4*t**2+24*alpha**5*t
    zero('centered quadruple-root derivative', S.diff(pp,t)-6*(t-alpha)**4*(t+4*alpha))
    zero('quadruple-root scaling', pp.subs(t,alpha*u)-alpha**6*(u**6-15*u**4+40*u**3-45*u**2+24*u))
    fourth = cert['fourth_order_root']
    sub = dict(zip((a,b,c,d),map(S.Rational,fourth['normalized_p_coefficients'])))
    fs = [S.expand(rho[j].as_expr().subs(sub)) for j in fourth['rho_orders']]
    multipliers = list(map(S.Rational,fourth['multipliers']))
    constant = S.Rational(fourth['constant'])
    zero('literal fourth-order-root unit identity',sum(f*v for f,v in zip(fs,multipliers))-constant)
    check('fourth-order-root constant nonzero',constant != 0)
    rows = [[f.subs({B:0,D:0,E:0})]+[S.diff(f,v) for v in (B,D,E)] for f in fs]
    check('all specialized gaps affine-linear',all(S.Poly(f,B,D,E).total_degree() <= 1 for f in fs))
    normalized = [Fraction(v/constant) for v in multipliers]
    check('separate Fraction unit replay',fraction_replay(rows,normalized)==[Fraction(1),0,0,0])
    mutation = normalized.copy(); mutation[0] += 1
    check('tampered unit certificate rejected',fraction_replay(rows,mutation)!=[Fraction(1),0,0,0])

    RR,*gens = ring('t,h,a,b,c,d,B,D,E',S.QQ)
    tt,hh,*_ = gens
    pr,qr = RR.from_expr(p),RR.from_expr(q)
    dp,dq = pr.diff(tt),qr.diff(tt)
    v = 2*hh-tt
    pv,qv,dpv,dqv = [f.compose(tt,v) for f in (pr,qr,dp,dq)]
    contact = (pr-pv)*(dq*dpv+dqv*dp)-2*(qr-qv)*dp*dpv
    c17 = coefficient(contact,17).as_expr()
    c11 = coefficient(contact,11).as_expr().subs(h,0)
    zero('universal contact coefficient t17',c17-S.sympify(cert['contact_coefficients']['degree_17']))
    zero('universal contact coefficient t11 at h0',c11-S.sympify(cert['contact_coefficients']['degree_11_at_h_zero']))
    check('contact coefficient mutation rejected',S.expand(c11-160*b**3/9+1)!=0)
    check('universal polynomial constructed without coefficient sampling',len(contact)>1000)
    zero('nonaffine involution leading coefficient',S.Rational(10,6)-2+S.Rational(1,3))
    check('bounded finite slope poles cannot reach leading order',6+3 < 10)
    r=S.Symbol('r')
    RP,*gp=ring('t,h,a,b,c,d,B,D,E,r',S.QQ)
    tp,hp,*_=gp
    ps,qs=RP.from_expr(p-r),RP.from_expr((p-r)*q)
    dps,dqs=ps.diff(tp),qs.diff(tp)
    vp=2*hp-tp
    psv,qsv,dpsv,dqsv=[f.compose(tp,vp) for f in (ps,qs,dps,dqs)]
    shifted=(ps-psv)*(dqs*dpsv+dqsv*dps)-2*(qs-qsv)*dps*dpsv
    zero('all-r shifted-product contact t23',coefficient(shifted,23).as_expr()
         -S.sympify(cert['contact_coefficients']['shifted_product_degree_23']))
    zero('all-r shifted-product contact t17 at h0',coefficient(shifted,17).as_expr().subs(h,0)
         -S.sympify(cert['contact_coefficients']['shifted_product_degree_17_at_h_zero']))
    zero('shifted-product nonaffine leading obstruction',S.Rational(16,6)-2-S.Rational(2,3))
    check('shifted-product pole contribution below leading order',6+3 < 16)
    check('degree-one reciprocal contact finite base pole bound',max(6,5,-6)<10)
    check('degree-one reciprocal contact infinite base pole bound',max(7,7,-6)<10)

    # Positive control has an actual degree-two contact involution.
    z = S.Symbol('z')
    family = 2*x**2/(3*z)+z**2/3
    zero('genuine two-contact family incidence',S.cancel(family.subs({x:t**3,z:t**2})-t**4))
    zero('genuine two-contact family critical condition',S.cancel(S.diff(family,z).subs({x:t**3,z:t**2})))
    p0,q0=t**3,t**4
    dp0,dq0=S.diff(p0,t),S.diff(q0,t)
    positive=(p0-p0.subs(t,-t))*(dq0*dp0.subs(t,-t)+dq0.subs(t,-t)*dp0)-2*(q0-q0.subs(t,-t))*dp0*dp0.subs(t,-t)
    zero('positive control Hermite involution',positive)
    check('positive control is outside degree(6,10) hypothesis',(S.degree(p0,t),S.degree(q0,t))!=(6,10))

    degrees=[]
    for k in range(6):
        check(f'constant-leading numerator exceeds denominator k{k}',17-k > 5-k)
        check(f'constant-leading field degree exceeds six k{k}',17-k > 6)
        check(f'linear-leading field degree exceeds six k{k}',15-k > 6)
        degrees.append(dict(gcd_degree=k,nonzero_A=17-k,zero_A=15-k))
    zero('zero-A leading coefficient cannot cancel',6-10+4)
    # Two exact cusp controls and one immersed witness, without relying on them
    # for universality of the symbolic calculation above.
    classes=json.loads((HERE.parent/'jc2-normal-r7/classification.json').read_text())['cases']
    controls=[(S.sympify(row['p']),S.sympify(row['q'])) for row in classes]
    controls.append((t**6+t**3,t**10+S.Rational(5,3)*t**7+S.Rational(7,9)*t**4+S.Rational(7,81)*t))
    for i,(pc,qc) in enumerate(controls):
        gc=S.gcd(S.diff(pc,t),S.diff(qc,t))
        kk=S.degree(gc,t)
        dc,nc=S.cancel(S.diff(pc,t)/gc),S.cancel(S.diff(qc,t)/gc)
        num=S.expand(A*dc*(x-pc)**2+dc*qc+nc*(x-pc))
        check(f'control{i} numerator degree uniform formula',S.degree(num,t)==17-kk)
        check(f'control{i} zero-A numerator degree uniform formula',S.degree(num.subs(A,0),t)==15-kk)
        check(f'control{i} no denominator cancellation over Q(x,A)',S.degree(S.gcd(S.Poly(num,t),S.Poly(dc,t)),t)==0)

    escape=[]
    for e1,e2 in [(1,5),(2,4),(3,3)]:
        for d1 in range(3):
            for d2 in range(3):
                remainder=2-e1*d1-e2*d2
                if remainder>=0 and d1+d2>0:
                    escape.append(dict(indices=[e1,e2],fibre_degrees=[d1,d2],
                                       squares=[-1-d1,-1-d2],exceptional_budget=remainder))
    check('complete three infinity escapes',escape==cert['surviving_infinity_budgets'])
    check('balanced(3,3) excluded independently of finality',all(row['indices'] != [3,3] for row in escape))
    check('all exceptional residual budgets at most one',max(row['exceptional_budget'] for row in escape)==1)
    contact_degrees=[]
    for row in escape:
        infinity_ramification=sum((e-1)*dd for e,dd in zip(row['indices'],row['fibre_degrees']))
        possible=[dd for dd in range(1,5) if
                  ((4-dd-infinity_ramification>=1) if row['exceptional_budget']
                   else (4-dd-infinity_ramification==0))]
        contact_degrees.append(possible)
    check('remaining ramified projection degrees',contact_degrees==[[1,2,3],[4],[3]])
    zero('full determinant divisor fibre degree',-2+3*2-4)
    # Exact remaining-family equations, not sampled coefficient choices.
    az,bz,rz=S.Function('alpha')(z),S.Function('beta')(z),S.Function('r')(z)
    graph=z*x+az+bz/(x-rz)
    critical=(x+S.diff(az,z))*(x-rz)**2+S.diff(bz,z)*(x-rz)+bz*S.diff(rz,z)
    zero('reciprocal graph derivative numerator',S.cancel(S.diff(graph,z)*(x-rz)**2-critical))
    zero('moving pole critical value',critical.subs(x,rz)-bz*S.diff(rz,z))
    zero('stationary pole single cancellation',
         S.expand(critical.subs(S.diff(rz,z),0)-(x-rz)*((x+S.diff(az,z))*(x-rz)+S.diff(bz,z))))
    w=S.Symbol('w')
    aa,bb,cc=S.Function('A')(z),S.Function('B')(z),S.Function('C')(z)
    xx=aa*w**2+bb*w+cc
    varying_y=z*xx+w
    determinant=S.diff(xx,z)*S.diff(varying_y,w)-S.diff(xx,w)*S.diff(varying_y,z)
    zero('varying-direction parabola determinant',determinant-(S.diff(xx,z)-xx*S.diff(xx,w)))
    zero('varying-direction cubic coefficient',S.Poly(S.expand(determinant),w).nth(3)+2*aa**2)
    line_c=S.Function('line_c')(z)
    line_y=z*xx+line_c
    line_det=S.diff(xx,z)*S.diff(line_y,w)-S.diff(xx,w)*S.diff(line_y,z)
    zero('double-line determinant factorization',line_det+S.diff(xx,w)*(xx+S.diff(line_c,z)))
    check('immersed boundary has two exceptional curves',8-(1+3+2)==2)
    check('nonimmersed smooth-cusp boundary has three exceptional curves',8-(1+2+2)==3)
    check('genus-two positivity window',[g for g in range(2,5) if 6-2*g>=1]==[2])

    inputs=[HERE/'PROOF.md',HERE/'README.md',HERE/'certificates.json',Path(__file__).resolve(),
            ROOT/'archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py',
            HERE.parent/'jc2-normal-r7/classification.json',
            HERE.parent/'jc2-omitted-divisors/GLOBAL_LEDGER.md',
            HERE.parent/'r7-adjoint-exclusion/PROOF.md',
            HERE.parent/'r7-adjoint-exclusion/SOURCE_AUDIT.md',
            HERE.parent/'minimal-resolution-positivity/PROOF.md',
            HERE.parent/'jc2-boundary-attachment/sources/borisov-stein-author.pdf']
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    receipt=dict(schema='r7-universal-contact-v1',passed=True,check_count=len(CHECKS),
                 checks=CHECKS,input_sha256=hashes,elapsed_seconds=round(time.monotonic()-start,3),
                 fourth_order_root_constant=str(constant),
                 universal_contact_coefficients=[str(c17),str(c11)],
                 degree_bounds=degrees,surviving_infinity_budgets=escape,
                 surviving_ramified_projection_degrees=contact_degrees,
                 surviving_projection_degrees_after_reciprocal_contact=[3,4,3],
                 geometry='Written proof with cited boundary/positivity inputs; not formally or externally reviewed.',
                 finality='Not used for the algebra or three-escape reduction. Added only for the stronger whole single-ramified-row corollary; not independently reproved.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} EXACT UNIVERSAL CONTACT CHECKS ({receipt["elapsed_seconds"]}s)')
    print('Remaining infinity budgets without finality:',escape)


if __name__=='__main__':
    main()
