#!/usr/bin/env python3
"""Prove that302 has a singleton rational isogeny class.

Three exact finite-field point enumerations and Mazur's prime-degree list.
No division-polynomial factorization, height search or rank upper bound.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(Path(__file__).resolve().parent))
from icarm_curve302 import GENERAL_WEIERSTRASS_COEFFICIENTS
PUBLIC=Path(__file__).with_name('icarm_curve302.py')
OUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_rational_isogeny_class_v1.json'
MAZUR_PRIMES=[2,3,5,7,11,13,17,19,37,43,67,163]


def prime(n):return n>1 and all(n%d for d in range(2,int(n**0.5)+1))


def build():
    assert all(a.denominator==1 for a in GENERAL_WEIERSTRASS_COEFFICIENTS)
    original=list(map(int,GENERAL_WEIERSTRASS_COEFFICIENTS));fibres=[]
    for p in [17,31,47]:
        assert prime(p);a1,a2,a3,a4,a6=[a%p for a in original];points=[]
        for x in range(p):
            for y in range(p):
                if (y*y+a1*x*y+a3*y-x*x*x-a2*x*x-a4*x-a6)%p:continue
                assert (a1*y-3*x*x-2*a2*x-a4)%p or (2*y+a1*x+a3)%p
                points.append([x,y])
        order=1+len(points);trace=p+1-order
        fibres.append({'prime':p,'reduced_general_ainvs':[a1,a2,a3,a4,a6],'affine_points':points,'order':order,'trace':trace,'frobenius_discriminant':trace*trace-4*p})
    assert [f['order'] for f in fibres]==[26,43,56]
    witnesses=[]
    for ell in MAZUR_PRIMES:
        assert prime(ell)
        for f in fibres:
            p=f['prime'];a=f['trace']
            if p==ell:continue
            values=[(x*x-a*x+p)%ell for x in range(ell)]
            if all(values):
                witnesses.append({'isogeny_prime':ell,'good_reduction_prime':p,'trace':a,
                                  'characteristic_polynomial_mod_ell_low_to_high':[p%ell,(-a)%ell,1],
                                  'all_polynomial_values':values});break
        else:raise AssertionError(f'No obstruction for degree {ell}')
    return {'schema':'curve302.rational-isogeny-class.v1','status':'PASS_SINGLETON_RATIONAL_ISOGENY_CLASS',
        'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),PUBLIC]},
        'theorem_source':{'author':'B. Mazur','title':'Rational Isogenies of Prime Degree','year':1978,'journal':'Inventiones mathematicae44,129–162','location':'Theorem1, pages129–130','url':'https://www.math.columbia.edu/~goldfeld/Mazur-Goldfeld1978.pdf','possible_prime_degrees':MAZUR_PRIMES},
        'finite_reductions':fibres,'prime_degree_obstructions':witnesses,
        'proof':'A rational cyclic ell-isogeny gives a Galois-stable line in E[ell]. At a good prime p different from ell, Frobenius must then have an eigenvalue in F_ell. Every possible ell from Mazur Theorem1 has a recorded irreducible polynomial X^2-a_p X+p, a contradiction. A rational isogeny of minimum degree from E to a nonisomorphic curve has cyclic kernel: otherwise its kernel contains E[n] for some n>1, allowing division by [n]. A nontrivial cyclic kernel has a rational prime-order subgroup, already excluded. Thus every Q-isogenous curve is Q-isomorphic to302.',
        'boundary':'This excludes rational isogenies to non-Q-isomorphic elliptic curves, not multiplication endomorphisms. It does not exclude isogenies over larger fields, recover a parent or prove an exact Mordell-Weil rank. All rational j-inverse family exclusions therefore remain valid if the proposed route is allowed a rational isogeny at specialization.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args();d=build()
    if args.check:assert d==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(d,sort_keys=True,indent=2,ensure_ascii=False)+'\n')
    print(d['status'])
