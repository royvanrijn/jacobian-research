#!/usr/bin/env python3
"""Exact finite CRT and genus transport checks for primitive x^2+5y^2.

Standard library only. This verifies finite algebra; it does not prove the
analytic finite-filter transfer lemma or any unrestricted global count.
"""


if not __debug__:
    raise RuntimeError(
        "This certificate checker requires active Python assertions; "
        "run without -O, -OO, or PYTHONOPTIMIZE."
    )
import argparse
from fractions import Fraction
from math import gcd
from pathlib import Path
import hashlib
import json

L = 40
TS = (1,2,5,10)
SIGNS = (1,-1)


def psi(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)


def chi5(n):
    return 0 if n % 5 == 0 else (1 if n % 5 in (1,4) else -1)


def chi(n):
    return psi(n)*chi5(n)


H = [a for a in range(L) if gcd(a,L) == 1 and chi(a) == 1]
A = {s:[a for a in H if psi(a) == s] for s in SIGNS}
assert A[1] == [1,9,21,29]
assert A[-1] == [3,7,23,27]

# First coefficient, second coefficient, first/second core genus signs,
# and the original target residues modulo 40 covered by this row.
ROWS = [
    (1,1,1,1,[2,10,18,22,30,38]),
    (1,5,1,1,[6,14,26,34]),
    (1,2,1,-1,[3,7,15,23,27,35]),
    (1,10,1,-1,[11,19,31,39]),
    (2,2,-1,-1,[12,20,28]),
    (2,10,-1,-1,[4,36]),
]


def count40(t,alpha,beta,target,u=1,v=1):
    return sum((u*a+t*v*b-target) % L == 0
               for a in A[alpha] for b in A[beta])


def crt_checks():
    all_counts = {}
    gamma_checks = 0
    for t in TS:
        period = t*L
        ht = [a for a in range(period)
              if gcd(a,period) == 1 and chi(a) == 1]
        assert len(ht) == t*len(H) == 8*t
        for alpha in SIGNS:
            first = [a for a in ht if psi(a) == alpha]
            assert len(first) == 4*t
            for beta in SIGNS:
                counts = [0]*period
                for a in first:
                    for b in A[beta]:
                        counts[(a+t*b) % period] += 1
                for target in range(period):
                    assert counts[target] == count40(t,alpha,beta,target)
                    gamma_checks += 1
                normalization = Fraction(period,len(ht)*len(H))
                assert normalization == Fraction(5,8)
                assert normalization*sum(counts)/period == Fraction(1,4)
                all_counts[t,alpha,beta] = counts
        assert Fraction(period,len(ht)*len(H))*(len(ht)*len(H))/period == 1

    # Check the whole principal-norm residue set by an independent direct
    # enumeration of coordinate pairs primitive at both primes dividing400.
    modulus = 400
    h400 = [a for a in range(modulus)
            if gcd(a,modulus) == 1 and chi(a) == 1]
    assert len(h400) == 80
    strata_residues = set()
    for coefficient,sign in [(1,1),(2,-1),(5,1),(10,-1)]:
        strata_residues |= {(coefficient*a) % modulus
                           for a in h400 if psi(a) == sign}
    direct_residues = {(x*x+5*y*y) % modulus
                       for x in range(modulus) for y in range(modulus)
                       if gcd(gcd(x,y),10) == 1}
    assert strata_residues == direct_residues
    sum_support = {(a+b) % modulus for a in strata_residues
                   for b in strata_residues}
    local_targets = {n for n in range(modulus) if n % 8 in {2,3,4,6,7}}
    assert sum_support == local_targets

    lookup = {}
    table = []
    for first_coeff,second_coeff,alpha,beta,residues in ROWS:
        g = gcd(first_coeff,second_coeff)
        assert first_coeff//g == 1
        t = second_coeff//g
        assert t in TS
        entries = []
        for residue in residues:
            assert residue not in lookup
            lookup[residue] = (first_coeff,second_coeff,alpha,beta,g,t)
            seen = set()
            for target in range(residue,modulus,L):
                assert target % g == 0
                reduced = target//g
                assert gcd(reduced,t) == 1
                count = all_counts[t,alpha,beta][reduced % (t*L)]
                assert count in (2,4)
                seen.add(count)
            assert len(seen) == 1
            count = seen.pop()
            entries.append({'target_mod40':residue,'pair_count':count,
                            'Gamma':str(Fraction(5,8)*count)})
        table.append({'first_coefficient':first_coeff,
                      'second_coefficient':second_coeff,
                      'first_genus':alpha,'second_genus':beta,
                      'common_coefficient':g,'reduced_t':t,
                      'entries':entries})
    assert set(lookup) == {n for n in range(L) if n % 8 in {2,3,4,6,7}}

    crt_cases = []
    for target in sorted(local_targets):
        r,s,alpha,beta,g,t = lookup[target % L]
        reduced = target//g
        period = t*L
        a,b = next((a,b) for a in range(period)
                   if gcd(a,period) == 1 and chi(a) == 1 and psi(a) == alpha
                   for b in A[beta] if (a+t*b-reduced) % period == 0)
        count = all_counts[t,alpha,beta][reduced % period]
        crt_cases.append({'target_mod400':target,'coefficients':[r,s],
                          'reduced_target_mod_tL':reduced % period,
                          'core_residue_witness':[a,b],
                          'pair_count':count,'Gamma':str(Fraction(5,8)*count)})
    assert len(crt_cases) == 250

    transport_cases = 0
    for t in TS:
        for alpha in SIGNS:
            for beta in SIGNS:
                for target in range(L):
                    for u in H:
                        for v in H:
                            direct = count40(t,alpha,beta,target,u,v)
                            transformed = count40(t,alpha*psi(u),beta*psi(v),target)
                            assert direct == transformed
                            transport_cases += 1
                    for common in H:
                        assert count40(t,alpha,beta,common*target) == count40(
                            t,alpha*psi(common),beta*psi(common),target)
                        transport_cases += 1
    # Deliberately challenge the incorrect unchanged-sign gcd rule.
    assert count40(1,1,1,18) == 2
    assert count40(1,1,1,6) == 0
    assert count40(1,-1,-1,6) == 2
    # A wrong genus in the N=4 row must lose its local main.
    assert count40(5,1,-1,2) == 0
    assert count40(5,-1,-1,2) == 4

    deep_five_cases = 0
    for exponent in range(101):
        for unit_part in [2,3,4,6,7]:
            target = 5**exponent*unit_part
            r,s,alpha,beta,g,t = lookup[target % L]
            reduced = target//g
            assert target % g == 0 and gcd(reduced,t) == 1
            assert all_counts[t,alpha,beta][reduced % (t*L)] > 0
            if exponent:
                assert r % 5 and s % 5
            deep_five_cases += 1

    return {
        'scope':'Finite CRT identities and genus transport, not an unrestricted analytic theorem.',
        'base_period':L,'uniform_original_target_period':modulus,
        'H':H,'A_plus':A[1],'A_minus':A[-1],
        'principal_norm_residue_count_mod400':len(strata_residues),
        'admissible_target_residue_count_mod400':len(local_targets),
        'finite_Gamma_checks':gamma_checks,
        'cofactor_and_gcd_transport_checks':transport_cases,
        'deep_five_challenges':deep_five_cases,
        'minimum_selected_Gamma':'5/4',
        'maximum_selected_Gamma':'5/2',
        'stratum_table':table,
        'all_admissible_CRT_cases':crt_cases,
        'negative_checks':{'unchanged_genus_gcd_rule_rejected':True,
                           'wrong_genus_in_target_four_rejected':True},
    }




def full_pair_density(t,u,v,alpha,beta):
    """Independent symmetric definition of FV1 at period tL."""
    period=t*L
    units=[a for a in range(period)
           if gcd(a,period)==1 and chi(a)==1]
    counts=[0]*period
    for a in units:
        if psi(a)!=alpha:continue
        for b in units:
            if psi(b)==beta:
                counts[(u*a+t*v*b)%period]+=1
    return [Fraction(period*x,len(units)**2) for x in counts]


def full_density_transport_checks():
    """Symmetric pair counts, coefficients and inverse target scaling."""
    checks=0
    cofactors=[1,3,7,9,21,29]
    assert all(gcd(n,400)==1 and chi(n)==1 for n in cofactors)
    rows=[]
    for t in TS:
        period=t*L
        base={(alpha,beta):full_pair_density(t,1,1,alpha,beta)
              for alpha in SIGNS for beta in SIGNS}
        for (alpha,beta),density in base.items():
            for target in range(period):
                q=sum((target-t*b)%L in H
                      and psi((target-t*b)%L)==alpha
                      for b in A[beta])
                assert density[target]==Fraction(5,8)*q
                assert density[target]==density[target%L]
                if density[target]:assert gcd(target,t)==1
                checks+=1
            for u in cofactors:
                for v in cofactors:
                    transformed=full_pair_density(
                        t,u,v,alpha*psi(u),beta*psi(v))
                    assert transformed==density
                    checks+=1
            for common in cofactors:
                transformed=base[alpha*psi(common),beta*psi(common)]
                for target in range(period):
                    reduced=(target*pow(common,-1,period))%period
                    assert density[target]==transformed[reduced]
                    checks+=1
    for t,alpha,beta in ((1,1,1),(5,1,1),(2,1,-1),(10,1,-1),
                         (1,-1,-1),(5,-1,-1)):
        density=full_pair_density(t,1,1,alpha,beta)
        grouped={}
        for target in range(L):
            if density[target]:
                grouped.setdefault(str(density[target]),[]).append(target)
        rows.append(dict(t=t,alpha=alpha,beta=beta,densities=grouped))
    assert checks==20736
    return dict(
        scope="Finite symmetric FV1 density and genus/cofactor transport only.",
        checks=checks,cofactors=cofactors,rows=rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",type=Path,
        default=Path(__file__).with_name("primitive_norm_genus_filter_checks.json"))
    args=parser.parse_args()
    result=dict(
        scope="Exact finite genus-filter, CRT and cofactor identities. "
              "No analytic estimate, unweighted sieve theorem, numerical "
              "cutoff or unrestricted exact exception claim is certified.",
        full_density_transport=full_density_transport_checks(),
        local_CRT=crt_checks(),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    text=json.dumps(result,indent=2)+"\n"
    args.output.write_text(text)
    print(json.dumps(dict(
        full_density_transport_checks=result["full_density_transport"]["checks"],
        CRT_Gamma_checks=result["local_CRT"]["finite_Gamma_checks"],
        CRT_transport_checks=result["local_CRT"]["cofactor_and_gcd_transport_checks"],
        admissible_target_classes_mod400=result["local_CRT"]["admissible_target_residue_count_mod400"],
        deep_five_challenges=result["local_CRT"]["deep_five_challenges"],
        output=str(args.output)),indent=2))


if __name__=="__main__":
    main()
