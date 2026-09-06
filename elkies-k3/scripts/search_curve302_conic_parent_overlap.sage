#!/usr/bin/env sage-python
"""Bounded exact modular screen of 2898 conics on the pinned 302 pencil.

Each conic passes through four of the nine basepoints and one of P9,...,P31.
The question is equality of its degree-two cover with the existing K3 cover.
A clean unequal branch divisor (or nonsquare relative constant) excludes
that conic only. Unusable reductions and surviving rows remain UNKNOWN.
"""
import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from itertools import combinations
import json
from pathlib import Path
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, matrix
from sage.env import SAGE_VERSION

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'artifacts/generated-results/elkies-k3-curve302-nine-direction-k3-v1.json'
PUBLIC=ROOT/'elliptic-curves/cas/icarm_curve302.py'
OUTPUT=ROOT/'artifacts/generated-results/elkies-k3-curve302-conic-parent-overlap-v1.json'
PRIMES=(101,103,107,109,127)

def digest(p): return sha256(p.read_bytes()).hexdigest()

def build():
    old=json.loads(BASE.read_text()); c=old['pointed_cubic_family']
    public=SourceFileLoader('conic302public',str(PUBLIC)).load_module()
    E=EllipticCurve(QQ,[QQ(str(a)) for a in public.GENERAL_WEIERSTRASS_COEFFICIENTS])
    points=[(36*QQ(str(x))+15,108*(2*QQ(str(y))+QQ(str(x))+1),QQ(1)) for x,y in public.POINTS]
    bases=[tuple(map(QQ,c['zero']))]+[tuple(map(QQ,r)) for r in c['constant_sections']]
    jobs=[(b,i) for b in combinations(range(9),4) for i in range(8,31)]
    pending=set(jobs); decisions={}; skipped=Counter(); counts=[]
    for prime in PRIMES:
        k=GF(prime); R=PolynomialRing(k,names=('X','Y','Z'));x,y,z=R.gens()
        if E.discriminant().valuation(prime):
            raise ArithmeticError('declared screen prime is bad for target')
        def poly(records):
            return sum(k(QQ(r['coefficient']))*x**r['exponents'][0]*y**r['exponents'][1]*z**r['exponents'][2] for r in records)
        try:
            f0,f1=poly(c['F0']),poly(c['F1'])
            bp=[tuple(k(v) for v in p) for p in bases]
        except (ValueError,ZeroDivisionError):
            skipped['source_denominator']+=len(pending); continue
        U=PolynomialRing(k,'u');u=U.gen()
        n=list(map(lambda a:k(QQ(a)),c['N_coefficients_low_to_high']));d=list(map(lambda a:k(QQ(a)),c['D_coefficients_low_to_high']))
        n += [k(0)]*(3-len(n));d += [k(0)]*(3-len(d))
        target_disc=(n[1]-u*d[1])**2-4*(n[2]-u*d[2])*(n[0]-u*d[0])
        if target_disc.degree()!=2 or target_disc.gcd(target_disc.derivative()).degree()!=0:
            skipped['target_branch_bad_reduction']+=len(pending);continue
        # Positive scalar and negative branch controls for the comparison.
        def compare(disc):
            scalar=disc[2]/target_disc[2]
            if disc!=scalar*target_disc: return 'branch_divisor_mismatch'
            if not scalar.is_square(): return 'relative_constant_nonsquare'
            return None
        assert compare(target_disc)==None and compare(4*target_disc)==None
        assert compare(target_disc+1)=='branch_divisor_mismatch'
        S=PolynomialRing(k,'s');s=S.gen()
        monomials=[x*x,x*y,y*y,x*z,y*z,z*z]
        before=len(pending)
        for subset,index in sorted(pending):
            try: p=tuple(k(a) for a in points[index])
            except (ValueError,ZeroDivisionError):
                skipped['point_denominator']+=1;continue
            selected=[bp[i] for i in subset]+[p]
            A=matrix(k,[[m(*q) for m in monomials] for q in selected])
            if A.rank()!=5:
                skipped['conic_conditions_not_rank5']+=1;continue
            conic=sum(a*m for a,m in zip(A.right_kernel().basis()[0],monomials))
            hessian=matrix(k,[[conic.derivative(v).derivative(w) for w in (x,y,z)] for v in (x,y,z)])
            if not hessian.det():
                skipped['singular_conic_reduction']+=1;continue
            P=selected[0]
            if not P[2]: skipped['basepoint_at_infinity']+=1;continue
            X0,Y0=P[0]/P[2],P[1]/P[2]
            linear=conic.derivative(x)(X0,Y0,1)+s*conic.derivative(y)(X0,Y0,1)
            quadratic=S(conic(1,s,0))
            param=[X0*quadratic-linear,Y0*quadratic-s*linear,quadratic]
            assert S(conic(*param))==0
            g0,g1=S(f0(*param)),S(f1(*param))
            common=g0.gcd(g1)
            if not common: skipped['zero_restriction']+=1;continue
            g0,g1=g0//common,g1//common
            if max(g0.degree(),g1.degree())!=2 or g0.gcd(g1).degree()!=0:
                skipped['residual_degree_not_two']+=1;continue
            aa,bb,cc=(U(g0[j])+u*U(g1[j]) for j in (2,1,0))
            disc=bb*bb-4*aa*cc
            if disc.degree()!=2 or disc.gcd(disc.derivative()).degree()!=0:
                skipped['conic_branch_bad_reduction']+=1;continue
            # The selected rational point supplies a rational preimage of u=0.
            assert disc(0).is_square()
            reason=compare(disc)
            if reason:
                decisions[(subset,index)]={'basepoint_indices_one_based':[i+1 for i in subset],
                  'public_point_index_one_based':index+1,'prime':prime,'reason':reason,
                  'conic_coefficients':list(map(int,[conic.monomial_coefficient(m) for m in monomials])),
                  'branch_coefficients_low_to_high':[int(disc[j]) for j in range(3)],
                  'target_branch_coefficients_low_to_high':[int(target_disc[j]) for j in range(3)]}
        pending -= set(decisions)
        counts.append({'prime':prime,'entered':before,'excluded':before-len(pending),'remaining':len(pending)})
        print('CONIC302|p=%s|excluded=%s|remaining=%s'%(prime,before-len(pending),len(pending)),flush=True)
        if not pending: break
    return {'schema':'elkies-k3.curve302-conic-parent-overlap.v1','status':'PASS_BOUNDED_MODULAR_SCREEN',
      'scope':{'basepoint_subsets':126,'public_points_one_based':[9,31],'conics':len(jobs),'prime_budget':list(PRIMES)},
      'excluded':len(decisions),'unresolved':len(pending),'prime_counts':counts,'skipped_reduction_counts':dict(sorted(skipped.items())),
      'decisions':[decisions[j] for j in sorted(decisions)],
      'unresolved_rows':[{'basepoints':[i+1 for i in b],'public_point':i+1} for b,i in sorted(pending)],
      'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in (BASE,PUBLIC)},
      'source_sha256':digest(Path(__file__)),'software':{'sage':SAGE_VERSION},
      'claim_boundary':'Only these conics are screened for the existing degree-two cover. No exclusion of other K3s, original parents, multisections of higher degree, point combinations, or generic MW rank is implied. Unresolved reductions are not exclusions.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    text=json.dumps(build(),indent=2,sort_keys=True)+'\n'
    if args.check: assert OUTPUT.read_text()==text
    else: OUTPUT.write_text(text)
if __name__=='__main__': main()
