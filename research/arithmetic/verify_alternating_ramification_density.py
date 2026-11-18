#!/usr/bin/env python3
"""Exact controls for zero-branch obstructions to the rank-256 Hasse family.

Run: .venv/bin/python research/arithmetic/verify_alternating_ramification_density.py
The arithmetic density and Hensel arguments are in the accompanying note.
"""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
from collections import Counter
from itertools import permutations,product
import sympy as sp
from verify_rational_degree_ten_pencil import rabin_irreducible
from verify_alternating_subgroup_hasse import permutation,hits,generated_group,orbit_masks

def cycle_type(value):
    unseen=set(range(len(value)));out=[]
    while unseen:
        start=min(unseen);point=start;length=0
        while point in unseen:
            unseen.remove(point);length+=1;point=value[point]
        out.append(length)
    return tuple(sorted(out))

def span(vectors):
    out={0}
    for value in vectors:
        out|={value^old for old in tuple(out)}
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    assert sp.__version__=="1.14.0"
    z,u,s=sp.symbols("z u s")
    H=9*z**10-864*z**8+28161*z**6-350150*z**4+1127500*z**2
    R=9*u**3-414*u**2+1836*u+400
    g=R.subs(u,z**2)
    assert sp.expand(H+250000-(z*z-25)**2*g)==0
    cubic_disc=sp.discriminant(R,u)
    sextic_disc=sp.discriminant(g,z)
    assert cubic_disc==2**6*3**6*5**5*13**2*17
    assert sextic_disc==-2**22*3**14*5**12*13**4*17**2
    assert g.subs(z,5)==g.subs(z,-5)==-5**2*13**2*17
    cubic7=sp.Poly(R,u,modulus=7).monic()
    assert cubic7==sp.Poly(u**3+3*u**2+u+4,u,modulus=7)
    assert rabin_irreducible(cubic7)
    assert rabin_irreducible(sp.Poly(g,z,modulus=7))
    q2=sp.Poly(z*z+3,z,modulus=29)
    q4=sp.Poly(z**4+9*z*z+3,z,modulus=29)
    assert sp.Poly(g,z,modulus=29).monic()==q2*q4
    assert rabin_irreducible(q2) and rabin_irreducible(q4)
    assert sp.gcd(q2,q4).degree()==0
    assert -sp.Rational(400,9)==-sp.Rational(20,3)**2
    assert -sp.Rational(-414,9)==46
    assert cubic_disc>0
    assert span([0b001,0b010,0b100])==set(range(8))
    assert span([0b110,0b101,0b011])=={0,3,5,6}
    for odd in [1,2,4,7]:
        assert span([0b110,0b101,odd])==set(range(8))

    wreath=[]
    for sigma in permutations(range(3)):
        for bits in product((0,1),repeat=3):
            value=tuple(2*sigma[j//2]+((j%2)^bits[j//2]) for j in range(6))
            wreath.append(value)
    assert len(wreath)==len(set(wreath))==48
    counts=Counter(cycle_type(value) for value in wreath)
    even_counts={",".join(map(str,t)):count for t,count in counts.items()
                 if all(length%2==0 for length in t)}
    assert even_counts=={"2,2,2":7,"2,4":6,"6":8}
    assert sum(even_counts.values())==21
    assert sp.Rational(21,48)==sp.Rational(7,16)

    d=permutation((0,1),(2,3))
    local=[]
    for cycles,expected in [
        ([(4,5,6,7,8,9)],12),
        ([(4,5),(6,7,8,9)],8),
        ([(4,5),(6,7),(8,9)],4)
    ]:
        residual=permutation(*cycles)
        odd=sum(len(c)-1 for c in cycles)%2
        h=permutation(*([(0,1)] if odd else []),*cycles)
        D=generated_group([d,h])
        assert len(D)==expected
        assert hits([d,h])==[]
        assert hits([d,h],True)==["H2"]
        orbits=sorted(mask.bit_count() for mask in orbit_masks([d,h]))
        assert all(length%2==0 for length in orbits)
        local.append(dict(residual_cycle_type=cycle_type(tuple(residual[j]-4 for j in range(4,10))),
            group_order=len(D),D_orbits=orbits,inertia=list(d),frobenius=list(h)))

    # Three-point upper bound for completely killing order-two inertia by P1 base change.
    base_change_bounds=[]
    for degree in range(1,21):
        max_even_fibres=0 if degree%2 else max(0,(4*degree-4)//degree)
        assert max_even_fibres<=3
        base_change_bounds.append(dict(degree=degree,upper_bound_fully_even_fibres=max_even_fibres))

    report=dict(schema_version=1,result="PASS",kind="exact_algebra_and_finite_group_controls",
        command=".venv/bin/python research/arithmetic/verify_alternating_ramification_density.py",
        python=platform.python_version(),sympy=sp.__version__,
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        zero_branch_residual=str(g),cubic_discriminant=int(cubic_disc),sextic_discriminant=int(sextic_disc),
        excluded_primes=[2,3,5,13,17],
        exact_modular_certificates={"7":"R/9=u^3+3u^2+u+4 irreducible; g/9 irreducible of degree six",
            "29":"g/9=(z^2+3)(z^4+9z^2+3), both factors irreducible"},
        full_signed_pair_group_order=48,
        signed_pair_cycle_type_counts={",".join(map(str,k)):v for k,v in sorted(counts.items())},
        all_even_cycle_counts=even_counts,bad_prime_density="7/16",
        zero_branch_local_groups=local,base_change_degree_sanity=base_change_bounds,
        negative_controls=["All individual elements are covered by H1,H3,W5 but each displayed D is missed.",
            "Adding H2 repairs all three displayed decomposition groups.",
            "Finite congruence restrictions cannot control the infinite positive-density prime set."],
        proof_scope="Finite checker proves algebra, modular irreducibility and group counts. Full wreath Galois group also uses the written real-root and complex-conjugation argument. Infinite bad primes, integer density zero and progression obstructions import Chebotarev and integral HIT, with the written Hensel/CRT proof. The all-degree base-change bound is Riemann-Hurwitz, not a consequence of its finite sanity table.",
        independent_replay=False,formal_verification=False,external_review=False)
    output=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(output)
    print(output,end="")

if __name__=="__main__":
    main()
