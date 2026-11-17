"""Integer-only necessary filter for actual affine/Keller omitted primes.

A valid signed key sequence determines Laurent key forms, not a realization
as an affine boundary or Keller map. Nonnegative beta_0 at every recurrence
is excluded by the polynomial-chain theorem; a negative beta_0 only passes
this necessary filter. No key polynomial is expanded here.
"""
from math import gcd


def analyze(weights):
    w = tuple(weights)
    if len(w) < 2 or any(type(a) is not int for a in w) or w[0] <= 0:
        raise ValueError("integer weights, length >=2 and omega_0>0 required")
    d = [w[0]]
    for a in w[1:]:
        d.append(gcd(d[-1], a))
    if d[-1] != 1:
        raise ValueError("normalized divisorial key sequence must have final gcd 1")
    alpha = [None] + [d[i-1]//d[i] for i in range(1, len(w))]
    rows = []
    for j in range(1, len(w)-1):
        if w[j+1] >= alpha[j]*w[j]:
            raise ValueError("strict key-weight descent violated at index " + str(j))
        remaining = alpha[j]*w[j]
        beta = [0]*j
        for i in range(j-1, 0, -1):
            assert remaining % d[i] == 0
            if alpha[i] != 1:
                beta[i] = ((remaining//d[i]) *
                           pow(w[i]//d[i], -1, alpha[i])) % alpha[i]
            remaining -= beta[i]*w[i]
        assert remaining % w[0] == 0
        beta[0] = remaining//w[0]
        assert sum(b*a for b,a in zip(beta,w)) == alpha[j]*w[j]
        rows.append(dict(index=j, next_key_index=j+1, alpha=alpha[j], beta=beta))
    escape = next((r for r in rows if r['beta'][0] < 0), None)
    polynomial = escape is None
    case = None if not polynomial else (
        "constants" if w[-1] > 0 else
        "one_generator" if w[-1] == 0 else "constant_mod_terminal")
    return dict(weights=list(w), gcds=d, alpha=alpha[1:], recurrences=rows,
                all_key_forms_polynomial=polynomial,
                first_Laurent_key_index=None if escape is None else escape['next_key_index'],
                regular_algebra_case=case,
                excluded_as_actual_affine_or_Keller_omitted_prime=polynomial,
                affine_or_Keller_realization_certified=False)


if __name__ == '__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('weights', nargs='+', type=int)
    print(json.dumps(analyze(parser.parse_args().weights), indent=2))
