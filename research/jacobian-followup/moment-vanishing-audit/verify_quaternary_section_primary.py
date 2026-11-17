#!/usr/bin/env python3
"""Exact monomial/Betti/Hilbert inputs for the written section-primary proof."""
from pathlib import Path
from collections import Counter, deque
from functools import lru_cache
import hashlib
import json
import re
import time
import sympy as S

P = Path(__file__).resolve().parent
EXPECTED = {
    'small_dimension_minimal_geometry_certificate.json': 'fe46cfd4cec299770eade67bafeed3f6aa4e5d9f74f456fb37a67d4ec67d6be0',
    'three_factor_associated_geometry_certificate.json': 'a6fd844bfdb11c64b64cfab60d68cbfbd8a12d62ec70d1dbf8cc59549c9a8f6c',
    'three_factor_minimal_primary_stdout.txt': 'fbeb16906d94c831bb16d67d44872587fbb51bd6c76bdf7e9d739b4a39549f08',
}
ACCEPTED = ['unmixed_monomial_depth', 'quaternary_modular_regularity',
            'quaternary_monomial_torsion']
PARTIAL = ['unmixed_modular_depth', 'unmixed_modular_minimized_depth',
           'quaternary_modular_section', 'quaternary_modular_linear_saturation']


def sha(name):
    return hashlib.sha256((P / name).read_bytes()).hexdigest()


def clean(name):
    return '\n'.join(x.strip() for x in (P / name).read_text().splitlines()
                     if x.strip() and not x.lstrip().startswith('//'))


def numerator(text, start, end):
    v = list(map(int, re.findall(r'-?\d+', text.split(start)[1].split(end)[0])))
    assert v[-1] == 0
    return v[:-1]


def polynomial(v, t):
    return sum(a*t**i for i, a in enumerate(v))


def betti(text, end):
    rows = [list(map(int, line.split())) for line in text.split('BETTI\n')[1].split(end)[0].splitlines()]
    assert rows and all(len(row) == len(rows[0]) for row in rows)
    assert rows[0][0] == 1 and all(a >= 0 for row in rows for a in row)
    return rows


def betti_numerator(rows, t):
    return S.expand(sum(a*(-1)**i*t**(r+i) for r, row in enumerate(rows)
                        for i, a in enumerate(row)))


def divides(a, b):
    return all(x <= y for x, y in zip(a, b))


def minimal(gens):
    out = []
    for m in sorted(set(gens), key=lambda v: (sum(v), v)):
        if not any(divides(g, m) for g in out):
            out.append(m)
    return out


def monomial(text):
    v = [0]*5
    pieces = re.findall(r'([abdef])(\d*)', text)
    assert ''.join(c+n for c, n in pieces) == text, text
    for c, n in pieces:
        v['abdef'.index(c)] += int(n or 1)
    return tuple(v)


def main():
    start = time.monotonic()
    for name, expected in EXPECTED.items():
        assert sha(name) == expected, name
    for stem in ACCEPTED + PARTIAL:
        rec = json.loads((P / (stem+'_run.json')).read_text())
        assert rec['terminal'] and rec['input_sha256'] == sha(stem+'.sing')
        assert rec['log_sha256'] == sha(stem+'_stdout.txt')
        if stem in ACCEPTED:
            assert rec['completed'] and not rec['timed_out'] and rec['exit_code'] == 0
            assert '?' not in (P / (stem+'_stdout.txt')).read_text()
        else:
            assert not rec['completed'] and rec['timed_out']
    t = S.symbols('t')
    old = json.loads((P / 'small_dimension_minimal_geometry_certificate.json').read_text())
    nh = polynomial(old['quaternary_generic_null_section']['first_Hilbert_numerator_over_5_variables'], t)
    oldlog = clean('three_factor_minimal_primary_stdout.txt')
    nu = polynomial(numerator(oldlog, 'UNMIXED_HILBERT', 'END'), t)
    # The characteristic-zero monomial resolution is bound to the retained U.
    src = (P / 'unmixed_monomial_depth.sing').read_text()
    utext = oldlog.split('UNMIXED_BASIS\n')[1].split('FULL_IN_UNMIXED')[0].strip()
    assert src.split('ideal U=')[1].split(';')[0].strip() == utext
    assert 'ring r=0,' in src and 'mres(L,0)' in src
    ub = betti(clean('unmixed_monomial_depth_stdout.txt'), 'END')
    assert betti_numerator(ub, t) == nu
    ureg = max(i for i, row in enumerate(ub) if any(row))
    assert ureg == 34
    # No depth>=2 conclusion follows: the initial quotient has projective dimension6.
    assert len(ub[0])-1 == 6 and any(row[-1] for row in ub)
    finite = S.cancel((nh-nu)/(1-t)**5)
    assert finite.is_polynomial(t)
    finite = S.Poly(finite, t)
    assert finite.eval(1) == 1635 and finite.degree() == 29
    assert all(a >= 0 for a in finite.all_coeffs())
    reg_bound = max(ureg, finite.degree())
    assert reg_bound == 34

    # The full modular run timed out only after a closed standard-basis phase.
    # Check its source and all checkpoint boundaries before reusing that phase.
    raw = (P / 'quaternary_modular_section_stdout.txt').read_text()
    assert '?' not in raw and 'SECTION_HILBERT' in raw and 'SATURATION_START' in raw
    assert 'SATURATION_HILBERT' not in raw
    original = (P / 'ternary_quaternary_null_section.sing').read_text().split('print("SECTION_STANDARD_START")')[0]
    modular = (P / 'quaternary_modular_section.sing').read_text()
    expected_prefix = original.replace('ring r=0,', 'ring r=32003,').replace('ring r1=0,', 'ring r1=32003,')
    assert modular.startswith(expected_prefix)
    assert 'ideal G=std(H);' in modular
    assert polynomial(numerator(raw, 'SECTION_HILBERT', 'SECTION_BASIS'), t) == nh
    basis = raw.split('SECTION_BASIS\n')[1].split('SATURATION_START')[0].strip()
    for stem in ['quaternary_modular_regularity', 'quaternary_monomial_torsion']:
        assert (P / (stem+'.sing')).read_text().split('ideal G=')[1].split(';')[0].strip() == basis
    rb = betti(clean('quaternary_modular_regularity_stdout.txt'), 'REGULARITY_BOUND')
    assert betti_numerator(rb, t) == nh
    assert max(i for i, row in enumerate(rb) if any(row)) == 38
    # Recompute the monomial origin saturation, rather than trusting its label.
    generators = []
    for line in basis.splitlines():
        m = re.match(r'[+-]?\d*([abdef]\d*(?:[abdef]\d*)*)', line)
        assert m, line[:100]
        generators.append(monomial(m.group(1)))
    L = minimal(generators)
    saturation = [(0,)*5]
    for i in range(5):
        piece = minimal([g[:i]+(0,)+g[i+1:] for g in L])
        saturation = minimal([tuple(max(x, y) for x, y in zip(a, b))
                              for a in saturation for b in piece])
    text = clean('quaternary_monomial_torsion_stdout.txt')
    reported = [monomial(x.rstrip(',')) for x in text.split('LEADING_SATURATION_BASIS\n')[1].split('END')[0].splitlines()]
    assert set(saturation) == set(reported)
    assert polynomial(numerator(text, 'LEADING_HILBERT', 'SATURATION_HILBERT'), t) == nh
    nsat = polynomial(numerator(text, 'SATURATION_HILBERT', 'LEADING_SATURATION_BASIS'), t)
    ht = S.Poly(S.cancel((nh-nsat)/(1-t)**5), t)
    @lru_cache(None)
    def in_L(m):
        return any(divides(g, m) for g in L)
    todo = deque(g for g in saturation if not in_L(g)); seen = set(todo)
    while todo:
        m = todo.popleft()
        for i in range(5):
            mm = m[:i]+(m[i]+1,)+m[i+1:]
            if mm not in seen and not in_L(mm):
                seen.add(mm); todo.append(mm)
        assert len(seen) <= 3000
    counts = Counter(map(sum, seen))
    assert len(seen) == 2582 and max(counts) == 37
    assert polynomial([counts[i] for i in range(38)], t) == ht.as_expr()
    upper = sum(a for d, a in counts.items() if d <= reg_bound)
    assert upper == 2576
    # Each selected isolator avoids its chosen prime and meets all the others.
    a,b,d,e,f = S.symbols('a b d e f'); variables = [a,b,d,e,f]
    c = -a-b-2*d-2*e-2*f
    primes = [[a,d,e],[b,d,f],[c,e,f],[b,c,f],[a,c,e],[a,b,d]]
    isolators = [b*c,a*c,a*b,d*e,d*f,e*f]
    for i, pol in enumerate(isolators):
        for j, prime in enumerate(primes):
            rem = S.groebner(prime, *variables).reduce(pol)[1]
            assert (rem != 0) == (i == j)
    names = list(EXPECTED)+['ternary_quaternary_null_section.sing',
             'quaternary_section_source_audit.json']
    for stem in ACCEPTED+PARTIAL:
        names += ['run_'+stem+'.py',stem+'.sing',stem+'_stdout.txt',stem+'_run.json']
    cert = {
        'claim': 'G4SECTIONPRIMARY20261114', 'field': 'C, characteristic zero',
        'section': 'I10 modulo c=-a-b-2d-2e-2f, in C[a,b,d,e,f]',
        'regularity_bound': reg_bound, 'origin_torsion_top_degree_bound': reg_bound,
        'primary_components': ['H:(bc)^infinity','H:(ac)^infinity','H:(ab)^infinity',
            'H:(de)^infinity','H:(df)^infinity','H:(ef)^infinity',
            'H+(a^35,b^35,d^35,e^35,f^35)'],
        'primary_components_count': 7, 'origin_power': 35,
        'origin_torsion_length_interval': [1635, upper],
        'injected_defect_Hilbert_function': [int(finite.nth(i)) for i in range(30)],
        'monomial_torsion_Hilbert_function': [counts[i] for i in range(38)],
        'monomial_torsion_length': 2582,
        'monomial_torsion_top_degree': 37,
        'modular_section_prime': 32003,
        'modular_section_Hilbert_numerator': old['quaternary_generic_null_section']['first_Hilbert_numerator_over_5_variables'],
        'monomial_generators': [list(g) for g in L],
        'monomial_saturation_generators': [list(g) for g in saturation],
        'unmixed_initial_Betti_table_over_Q': ub,
        'section_initial_Betti_table_over_Fp': rb,
        'provenance': {name: sha(name) for name in names},
        'checkpoint_boundary': 'The modular full run is terminal incomplete at saturation. Only its earlier completed std(H), full printed basis and matching Hilbert numerator are reused; no saturation result from that run is asserted.',
        'proof_boundary': 'The seven-way intersection is proved by graded separation, not an expanded CAS intersection. Hilbert matching gives degreewise DVR flatness; finite matrix ranks and Groebner flat families lift upper bounds only. No modular primary component, exact torsion length or depth>=2 assertion is transferred.',
        'remaining': 'Exact generic null torsion; remaining N=4 exceptional strata and global primary/nilpotence; full N=3 embedded schemes; arbitrary cubic GMC and lower-power GVC.',
        'assurance': 'Internal written proof and exact CAS/symbolic/combinatorial inputs; no independent replay, external review or formal verification.'}
    path = P / 'quaternary_section_primary_certificate.json'
    path.write_text(json.dumps(cert, indent=2)+'\n')
    rec = {'claim': cert['claim'], 'completed': True,
           'primary_components': 7, 'origin_power': 35,
           'regularity_bound': 34, 'origin_torsion_length_interval': [1635,2576],
           'checker_sha256': sha(Path(__file__).name), 'certificate_sha256': sha(path.name),
           'elapsed_seconds': time.monotonic()-start, 'assurance': cert['assurance']}
    (P / 'quaternary_section_primary_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2))


if __name__ == '__main__':
    main()
