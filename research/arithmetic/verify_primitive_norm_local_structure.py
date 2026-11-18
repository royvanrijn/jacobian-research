#!/usr/bin/env python3
"""Finite exact checks for PRIMITIVE_NORM_LOCAL_STRUCTURE.md.

The unrestricted local results use the note's Hensel arguments. This checker
verifies finite residue formulas and bounded nonrepresentation certificates.
Run from the repository root with Python 3 (standard library only).
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json
import platform

HERE = Path(__file__).resolve().parent


def qvalue(form, x, y):
    a, b, c = form
    return a*x*x+b*x*y+c*y*y


def discriminant(form):
    a, b, c = form
    return b*b-4*a*c


def prime_symbol(a, p):
    z = pow(a % p, (p-1)//2, p)
    assert z != 0
    return 1 if z == 1 else -1


def binary_counts(form, p, modulus=None):
    modulus = p if modulus is None else modulus
    counts = Counter()
    for x in range(modulus):
        for y in range(modulus):
            if x % p or y % p:
                counts[qvalue(form, x, y) % modulus] += 1
    return counts


def residue_sumset(values, modulus):
    return {(a+b) % modulus for a in values for b in values}


def global_values(form, bound):
    # Only a=1 is used. Completing the square gives
    # |y|<=sqrt(4*bound/|D|), |x|<=sqrt(bound)+|b*y|/2.
    assert form[0] == 1
    limit_y = isqrt(4*bound//(-discriminant(form)))
    values = set()
    for y in range(-limit_y, limit_y+1):
        limit_x = isqrt(bound)+abs(form[1]*y)+2
        for x in range(-limit_x, limit_x+1):
            n = qvalue(form, x, y)
            if 0 < n <= bound and gcd(x, y) == 1:
                values.add(n)
    return values



def maximal_order_local_tables():
    """Full finite parameter checks; the all-depth proof remains mathematical."""
    def expect_rejection(action):
        try:
            action()
        except AssertionError:
            return
        raise AssertionError("corrupted maximal-order local assertion accepted")

    def dyadic_images(m, modulus):
        if m % 4 == 3:
            if (-m) % 8 == 1:
                return set(range(modulus)), set(range(modulus))
            odd = {n for n in range(modulus) if n % 2}
            return odd, {n for n in range(modulus) if n % 2 == 0}
        if m % 2:
            vals = {n for n in range(modulus)
                    if n % 4 == 1 or n % 8 == (m+1) % 8}
            sums = {n for n in range(modulus)
                    if n % 8 in {2,3,4,6,7}}
            return vals, sums
        k = m//2
        vals = {n for n in range(modulus)
                if n % 8 in {1,(1+2*k) % 8}
                or n % 16 in {(2*k) % 16,(2*k+4) % 16}}
        sums = ({n for n in range(modulus) if n % 16}
                if k % 4 == 1 else set(range(modulus)))
        return vals, sums

    dyadic_rows = []
    dyadic_sums = {}
    modulus = 128
    # Two ramified odd-m classes, four even-m classes, and both
    # unramified dyadic splitting types, all with complete target domains.
    for m in (1,5,2,6,10,14,3,7):
        form = ((1,1,(m+1)//4) if m % 4 == 3 else (1,0,m))
        counts = binary_counts(form,2,modulus)
        vals = set(counts)
        sums = residue_sumset(vals,modulus)
        expected_vals,expected_sums = dyadic_images(m,modulus)
        assert vals == expected_vals and sums == expected_sums, m
        dyadic_sums[m] = sums
        # The generic stabilization bound is at most 2^7 for all of
        # these maximal-order forms. Verify the exact next-level count.
        higher = binary_counts(form,2,2*modulus)
        roots = [sum(counts[a]*counts[(n-a)%modulus]
                     for a in range(modulus)) for n in range(modulus)]
        high_roots = [sum(higher[a]*higher[(n-a)%(2*modulus)]
                          for a in range(2*modulus))
                      for n in range(2*modulus)]
        assert all(high_roots[n] == 8*roots[n%modulus]
                   for n in range(2*modulus)), m
        dyadic_rows.append({
            'squarefree_m':m,'fundamental_D':discriminant(form),
            'form':list(form),'modulus':modulus,
            'primitive_values':sorted(vals),'selfsum_residue_count':len(sums),
            'next_level':2*modulus,'normalized_pair_densities_stabilize':True,
        })

    # Direct convolution for every negative fundamental discriminant
    # through -800 challenges the parameter reduction independently.
    fundamental = []
    for m in range(1,801):
        if any(m % (d*d) == 0 for d in range(2,isqrt(m)+1)):
            continue
        form = ((1,1,(m+1)//4) if m % 4 == 3 else (1,0,m))
        D = discriminant(form)
        if -D > 800:
            continue
        vals = set(binary_counts(form,2,modulus))
        expected_vals,expected_sums = dyadic_images(m,modulus)
        assert vals == expected_vals
        assert residue_sumset(vals,modulus) == expected_sums
        fundamental.append(D)
    assert len(fundamental) == 245

    odd_rows = []
    density_cases = 0
    sums_at_three = {}
    for p in (3,5,7,11):
        nonsquare = next(u for u in range(2,p) if prime_symbol(u,p) == -1)
        for kappa in (1,nonsquare):
            # These are local ramified models. They need not be globally
            # maximal at other primes, which are not tested in this row.
            form = (1,0,p*kappa)
            mod = p**3
            counts = binary_counts(form,p,mod)
            vals = set(counts)
            qr = {u*u % p for u in range(1,p)}
            expected_vals = {
                n for n in range(mod)
                if n % p in qr
                or (n % p == 0 and (n//p) % p
                    and ((n//p)*pow(kappa,-1,p)) % p in qr)
            }
            assert vals == expected_vals
            # Both valuation strata have exactly two unit square roots
            # per free coordinate: the actual binary count is 2*p^3.
            assert all(v == 2*mod for v in counts.values())
            sums = residue_sumset(vals,mod)
            if p == 3:
                expected_sums = {
                    n for n in range(mod)
                    if n % p or n % 9 == (6*kappa) % 9
                }
                sums_at_three[kappa] = sums
            elif p % 4 == 1:
                expected_sums = set(range(mod))
            else:
                expected_sums = {n for n in range(mod) if n % (p*p)}
            assert sums == expected_sums, (p,kappa)
            # Independent closed ramified pair-density normalization,
            # using only character counts in F_p; check every p^3 target.
            eps = prime_symbol(-1,p)
            density_table = []
            for n in range(mod):
                actual_pairs = sum((n-a) % mod in vals for a in vals)
                actual = Fraction(4*actual_pairs,mod)
                if n % p:
                    expected = (Fraction(1)-Fraction(eps,p)
                                -Fraction(2*(1+prime_symbol(n,p)),p*p))
                elif n % (p*p):
                    ch = prime_symbol((n//p)*pow(kappa,-1,p),p)
                    expected = ((Fraction(2*(p-1),p) if eps == 1 else Fraction(0))
                                +Fraction(p-2-eps-2*ch,p*p))
                else:
                    expected = Fraction(2*(p*p-1),p*p) if eps == 1 else Fraction(0)
                assert actual == expected, (p,kappa,n,actual,expected)
                assert (actual > 0) == (n in expected_sums)
                if n < p*p:
                    density_table.append(str(actual))
                density_cases += 1
            odd_rows.append({
                'prime':p,'kappa':kappa,'form':list(form),'modulus':mod,
                'primitive_value_count':len(vals),'selfsum_residue_count':len(sums),
                'binary_count_on_each_allowed_value':2*mod,
                'pair_density_by_target_mod_p2':density_table,
            })

    # These corruptions target precisely the delicate omitted cases.
    expect_rejection(lambda: assert_equal_local(0 in dyadic_sums[2],True))
    expect_rejection(lambda: assert_equal_local(0 in dyadic_sums[6],False))
    expect_rejection(lambda: assert_equal_local(3 in sums_at_three[1],True))
    expect_rejection(lambda: assert_equal_local(6 in sums_at_three[2],True))
    p5 = next(row for row in odd_rows if row['prime']==5 and row['kappa']==1)
    expect_rejection(lambda: assert_equal_local(p5['selfsum_residue_count'],125-5))
    return {
        'scope':'Finite parameter and normalization checks; full local images and all higher lifts are proved separately.',
        'dyadic_modulus':modulus,'dyadic_parameter_rows':dyadic_rows,
        'negative_fundamental_discriminants_checked':sorted(fundamental),
        'fundamental_discriminant_bound':800,'fundamental_discriminant_count':245,
        'odd_ramified_parameter_rows':odd_rows,
        'odd_ramified_pair_density_targets_checked':density_cases,
        'negative_controls':[
            'allowing16-divisible targets for k1 rejected',
            'excluding deep dyadic targets for k3 rejected',
            'omitting the ramified3 unit quotient rejected',
            'forgetting the kappa sign at3 rejected',
            'forbidding25-divisibility at ramified5 rejected',
        ],
    }


def assert_equal_local(actual,expected):
    assert actual == expected

def main():
    forms = [(1,0,1), (1,0,2), (1,0,5), (1,1,2),
             (1,1,3), (1,1,5), (2,1,3)]
    primes = [3,5,7,11,13,17,19,23,29]
    formula_cases = 0
    for p in primes:
        good = [f for f in forms if discriminant(f) % p]
        counts = {f: binary_counts(f, p) for f in good}
        for f in good:
            e1 = prime_symbol(discriminant(f), p)
            for g in good:
                e2 = prime_symbol(discriminant(g), p)
                for n in range(p):
                    actual = sum(counts[f][t]*counts[g][(n-t) % p]
                                 for t in range(p))
                    if n:
                        expected = p**3-(2+e1*e2)*p+e1+e2
                    else:
                        expected = (p**3+e1*e2*p**2
                                    -(2+e1+e2+e1*e2)*p+1+e1+e2)
                    assert actual == expected > 0, (p,f,g,n,actual,expected)
                    assert Fraction(actual, p**3) > 0
                    formula_cases += 1
    dyadic = []
    for modulus in [8,16,32,64,128]:
        vals = set(binary_counts((1,0,2), 2, modulus))
        sums = residue_sumset(vals, modulus)
        if modulus >= 16:
            assert sums == {n for n in range(modulus) if n % 16}
        dyadic.append({'D':-8, 'modulus':modulus,
                       'primitive_residues':sorted(vals),
                       'sumset_size':len(sums)})
        vals = set(binary_counts((1,0,5), 2, modulus))
        sums = residue_sumset(vals, modulus)
        assert sums == {n for n in range(modulus) if n % 8 in {2,3,4,6,7}}
        dyadic.append({'D':-20, 'modulus':modulus,
                       'primitive_residues':sorted(vals),
                       'sumset_size':len(sums)})
    assert set(binary_counts((1,0,2),2,16)) == {1,2,3,6,9,11}
    assert set(binary_counts((1,0,5),2,8)) == {1,5,6}
    ramified = []
    for p in [3,7,11,19]:
        modulus = p*p
        form = (1,1,(p+1)//4)
        vals = set(binary_counts(form,p,modulus))
        sums = residue_sumset(vals,modulus)
        if p == 3:
            expected = {n for n in range(modulus) if n % 3 or n % 9 == 6}
        else:
            expected = set(range(1,modulus))
        assert sums == expected, (p,sums ^ expected)
        ramified.append({'p':p,'modulus':modulus,'sumset_size':len(sums)})
    # Finite residue identities used by the additive-order corollary PL5.
    gaussian_r = set(binary_counts((1,0,1),2,8))
    assert gaussian_r == {1,2,5}
    gaussian_s = residue_sumset(gaussian_r,8)
    gaussian_three = {(a+b)%8 for a in gaussian_s for b in gaussian_r}
    assert gaussian_s == {2,3,4,6,7}
    assert gaussian_three == set(range(8))-{2}
    assert residue_sumset(gaussian_s,8) == set(range(8))
    eisenstein_r = set(binary_counts((1,1,1),3,9))
    assert eisenstein_r == {1,3,4,7}
    assert residue_sumset(eisenstein_r,9) == set(range(9))-{0,3}
    additive_residues = {
        'Gaussian_modulus':8,
        'Gaussian_primitive_norms':sorted(gaussian_r),
        'Gaussian_two_sums':sorted(gaussian_s),
        'Gaussian_three_sums':sorted(gaussian_three),
        'Gaussian_four_sums':sorted(residue_sumset(gaussian_s,8)),
        'Eisenstein_modulus':9,
        'Eisenstein_primitive_norms':sorted(eisenstein_r),
        'Eisenstein_two_sums':sorted(residue_sumset(eisenstein_r,9)),
    }
    # Exact PL2 stabilization at h0=2*v_p(D)+1, including a conductor prime.
    stabilized = []
    for form,p in [((1,0,2),2),((1,0,5),2),((1,1,1),3),((1,0,3),2)]:
        d = abs(discriminant(form))
        t = 0
        while d % p == 0:
            d //= p
            t += 1
        h0 = 2*t+1
        modulus = p**h0
        higher_modulus = modulus*p
        counts = binary_counts(form,p,modulus)
        higher_counts = binary_counts(form,p,higher_modulus)
        roots = [sum(counts[j]*counts[(n-j)%modulus]
                     for j in range(modulus)) for n in range(modulus)]
        higher_roots = [sum(higher_counts[j]*higher_counts[(n-j)%higher_modulus]
                            for j in range(higher_modulus))
                        for n in range(higher_modulus)]
        assert all(count == p**3*roots[n%modulus]
                   for n,count in enumerate(higher_roots))
        stabilized.append({'form':list(form),'prime':p,'h0':h0,
                           'next_level_target_residues':higher_modulus,
                           'normalized_counts_equal':True})
    # PL6 explicit bad factors: the -8 count at its proved stabilization level.
    modulus = 128
    c = binary_counts((1,0,2),2,modulus)
    beta2_table = []
    for n in range(modulus):
        actual = Fraction(sum(c[j]*c[(n-j)%modulus] for j in range(modulus)),
                          modulus**3)
        expected = (Fraction(0) if n % 16 == 0 else
                    Fraction(5,4) if n % 8 == 4 else Fraction(1,2))
        assert actual == expected
        if n < 16:
            beta2_table.append(str(actual))
    explicit_ramified_cases = 0
    for p in [3,7,11,19]:
        modulus = p*p
        c = binary_counts((1,1,(p+1)//4),p,modulus)
        for n in range(modulus):
            actual = Fraction(sum(c[j]*c[(n-j)%modulus] for j in range(modulus)),
                              modulus**3)
            if n % p:
                expected = Fraction(p*p+p-2*(1+prime_symbol(n,p)),p*p)
            elif n:
                expected = Fraction(p-1-2*prime_symbol(n//p,p),p*p)
            else:
                expected = Fraction(0)
            assert actual == expected, (p,n,actual,expected)
            explicit_ramified_cases += 1
        c2 = binary_counts((1,1,(p+1)//4),2)
        for n in range(2):
            actual = Fraction(sum(c2[j]*c2[(n-j)%2] for j in range(2)),8)
            expected = (Fraction(9,8) if n == 0 else Fraction(0)) if p % 8 == 3 else (
                        Fraction(5,8) if n == 0 else Fraction(1,2))
            assert actual == expected
    explicit_bad_factors = {'D_minus8_mod16':[str(x) for x in beta2_table],
                            'odd_ramified_factor_cases':explicit_ramified_cases}
    counterexamples = []
    for form,n in [((1,1,2),7),((1,1,5),4),((1,0,5),3),((1,0,5),4)]:
        vals = global_values(form,n-1)
        assert not any(n-a in vals for a in vals)
        counterexamples.append({'form':list(form),'N':n,
                                'all_positive_primitive_values_below_N':sorted(vals)})
    assert 2 not in global_values((1,1,6),2)
    # A conductor-prime regression: the maximal-order valuation<=1 rule fails.
    assert gcd(0,1) == 1 and qvalue((1,0,4),0,1) == 4
    assert gcd(2,1) == 1 and qvalue((1,0,4),2,1) == 8
    # Complete finite calculations for the exact higher-sumset table PL6.
    # PL5 proves the cutoff, so these bounded tables give exact exceptions.
    expected_three = {
        -3:[], -11:[], -19:[5], -43:[5,7,9,11,17,21],
        -67:[5,7,9,11,13,15,17,23,27,29,33,45],
        -163:[5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,
              47,51,53,57,59,61,65,67,69,71,75,77,79,81,93,111,117,121],
    }
    exact_higher = []
    for c in [1,3,5,11,17,41]:
        form = (1,1,c)
        d = discriminant(form)
        all_values = global_values(form,4*c-1)
        sums = {0}
        three_bad = []
        four_bad = []
        for k in range(1,5):
            sums = {a+b for a in sums for b in all_values if a+b < 4*c}
            if k == 3:
                three_bad = [n for n in range(3,3*c,2) if n not in sums]
                assert three_bad == expected_three[d], (d,three_bad)
            if k == 4:
                four_bad = [n for n in range(4,4*c,2) if n not in sums]
                assert four_bad == [n+1 for n in expected_three[d]], (d,four_bad)
        exact_higher.append({
            'discriminant':d,'c':c,'proved_three_cutoff':3*c,
            'all_primitive_norm_values_below_three_cutoff':
                sorted(n for n in all_values if n < 3*c),
            'three_summand_exceptions_at_correct_parity_above_minimum':three_bad,
            'proved_four_cutoff':4*c,
            'four_summand_exceptions_at_correct_parity_above_minimum':four_bad,
            'unrestricted_scope':'Exact by the elementary cutoff and translation proof PL5, not by an extrapolated scan.',
        })
    triple_regressions = []
    for form,d in [((1,0,1),-4),((1,1,2),-7),((1,0,2),-8)]:
        bound = 256
        vals = global_values(form,bound)
        three = {a+b for a in vals for b in vals if a+b <= bound}
        three = {a+b for a in three for b in vals if a+b <= bound}
        expected = {n for n in range(3,bound+1) if d != -4 or n % 8 != 2}
        assert three == expected, (d,sorted(three ^ expected))
        triple_regressions.append({'discriminant':d,'bound':bound,
                                    'matches_exact_three_sumset':True})
    conductor_residues = set(binary_counts((1,0,8),2,8))
    assert conductor_residues == {0,1,4}
    at_most_three = {0}
    for k in range(3):
        at_most_three |= {(a+b)%8 for a in at_most_three for b in conductor_residues}
    assert 7 not in at_most_three
    maximal_table = maximal_order_local_tables()
    receipt = {
        'scope':'Finite exact arithmetic checks; unrestricted local claims require the lifting proof.',
        'python':platform.python_version(),
        'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'good_prime_density_cases':formula_cases,
        'dyadic_tables':dyadic,
        'ramified_tables':ramified,
        'locally_admissible_global_counterexamples':counterexamples,
        'individual_norm_counterexample':{'D':-23,'n':2},
        'conductor_two_regressions':[4,8],
        'additive_order_residue_identities':additive_residues,
        'exact_density_stabilization_checks':stabilized,
        'exact_higher_sumsets_with_proved_cutoffs':exact_higher,
        'exact_triple_sumset_regressions':triple_regressions,
        'universal_order_four_sharpness_residues':sorted(conductor_residues),
        'explicit_bad_factors':explicit_bad_factors,
        'maximal_order_selfsum_parameter_checks':maximal_table,
        'outcome':'pass',
    }
    target = HERE/'primitive_norm_local_structure_checks.json'
    target.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:receipt[k] for k in ['outcome','good_prime_density_cases',
                      'locally_admissible_global_counterexamples']},indent=2))


if __name__ == '__main__':
    main()
