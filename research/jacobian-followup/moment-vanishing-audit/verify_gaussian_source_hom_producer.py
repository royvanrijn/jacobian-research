"""Replay retained source/Hom controls; --generate builds their CAS artifacts."""
from pathlib import Path
import argparse
import copy
import hashlib
import itertools
import json
import subprocess
import sys

import sympy as sp
import gaussian_source_hom_producer as G

p = Path(__file__).resolve().parent
proof = p/'WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md'
sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
dump = lambda x: json.dumps(x, indent=2)+'\n'
cases = [
    ('zero', dict(dimension=1, phase_terms=[]), True, 1, '0'),
    ('constant', dict(dimension=1, phase_terms=[[[0], '1']]), False, 1, None),
    ('quartic_whole_source', dict(dimension=1, phase_terms=[[[4], '1']]), False, 2, None),
    ('coupled_whole_source', dict(dimension=2, phase_terms=[[[1, 1], '1']]), False, 3, None),
    ('isotropic_cubic', dict(dimension=2, field=dict(kind='number_field', minimal_polynomial='a^2+1'),
                           phase_terms=[[[3, 0], '1'], [[2, 1], '3*a'], [[1, 2], '-3'], [[0, 3], '-a']]), True, 1, '1/3'),
    ('quartic_coefficient_field', dict(dimension=2, field=dict(kind='number_field', minimal_polynomial='a^4+1'),
                                     phase_terms=[[[1, 0], '1'], [[0, 1], 'a^2']]), True, 1, '1')]


def expectation(data, alpha, power):
    f, n, terms = G.normalize(data)
    xs = sp.symbols('x:'+str(n))
    phase = sum(f.expression(f.encode(c))*sp.prod(x**e for x, e in zip(xs, es)) for es, c in terms.items())
    expr = sp.Poly(sp.expand(sp.prod(x**e for x, e in zip(xs, alpha))*phase**power), *xs)
    out = f.K.zero
    for exponents, coefficient in expr.terms():
        out += f.parse(str(coefficient).replace('**', '^'))*f.K.convert(G.wick(exponents))
    return f, out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    args = parser.parse_args()
    controls, files = [], {}
    root = p/'gaussian_source_hom_examples';root.mkdir(exist_ok=True)
    actual_checks = 0
    results = {}
    for name, data, safe, rank, slope in cases:
        out = root/name
        if args.generate:
            out.mkdir(exist_ok=True)
            (out/'input.json').write_text(dump(data))
            result = G.decide(data, out/'cas', 20)
            (out/'result.json').write_text(dump(result))
        result = json.loads((out/'result.json').read_text())
        assert json.loads((out/'input.json').read_text()) == data
        assert G.verify_algebraic_parts(result)
        source = result['exact_unlocalized_source']
        f, n, phase = G.normalize(data)
        # Replay parsing of the retained CAS output, then independently replay
        # finite LEFT identities and the complete Hom bounds/linear system.
        parsed = G.parse_output((out/'cas/source.stdout').read_text(), f)
        G.check_integration_basis(parsed, n, f)
        G.certify_recurrence(parsed, f, n, phase)
        G.certify_coordinate_bounds(parsed, f)
        G.certify_original_ibp(parsed, n, phase, f)
        assert parsed == source
        for name_, hash_ in result['cas_files_sha256'].items():
            assert sha(out/'cas'/name_) == hash_
        assert result['all_original_marker_eventual_safety'] is safe
        assert source['rank'] == rank
        if safe:
            cutoff = result['all_marker_cutoff']
            assert cutoff['slope'] == slope and cutoff['all_positive_pure']
            C, E = sp.Rational(cutoff['slope']), sp.Rational(cutoff['intercept'])
            for alpha in itertools.product(range(6), repeat=n):
                if sum(alpha) > 5:
                    continue
                for m in range(5):
                    if m > C*sum(alpha)+E:
                        _, value = expectation(data, alpha, m)
                        assert not value
                        actual_checks += 1
        controls.append(dict(name=name, safe=safe, whole_module_generators=rank,
                             integration_integer_roots=source['integration_integer_roots'],
                             polynomial_hom_dimension=result['polynomial_hom']['dimension'],
                             all_marker_cutoff=result['all_marker_cutoff']))
        results[name] = result
        for f_ in out.rglob('*'):
            if f_.is_file():
                files[str(f_.relative_to(p))] = sha(f_)

    # Polynomial Hom controls with genuinely coupled generator values.
    D = [[[0, 1], '1']];one = [[[0, 0], '1']];minus_one = [[[0, 0], '-1']]
    Tminus = [[[1, 0], '-1']]
    synthetic = dict(rank=2, relation_columns=3,
                     relations=[[D, [], one], [minus_one, D, Tminus]],
                     coordinate_annihilators=[[[[0, 2], '1']], D],
                     coordinate_annihilator_left_witnesses=[[D, one, []], [[], one, []]])
    f = G.Field()
    G.certify_coordinate_bounds(synthetic, f)
    h = G.polynomial_hom(synthetic, f)
    assert h['dimension'] == 1 and h['generator_images'] == [[[[[1], '1']], [[[0], '1']]]]
    # Removing the third relation admits an extra independent constant map.
    larger = copy.deepcopy(synthetic);larger['relation_columns'] = 2
    larger['relations'] = [row[:2] for row in larger['relations']]
    larger['coordinate_annihilator_left_witnesses'] = [row[:2] for row in larger['coordinate_annihilator_left_witnesses']]
    G.certify_coordinate_bounds(larger, f)
    assert G.polynomial_hom(larger, f)['dimension'] == 2
    abstract_controls = ['coupled tuple (t,1)', 'missing cross-relation adds a polynomial map']
    for label, operator, dimension in [
        ('exact exponential source', [[[0, 1], '1'], [[0, 0], '-1']], 0),
        ('spurious larger exponential source', [[[0, 2], '1'], [[0, 1], '-1']], 1),
        ('origin torsion is retained', [[[1, 0], '1']], 0)]:
        s = dict(rank=1, relation_columns=1, relations=[[operator]], coordinate_annihilators=[operator],
                 coordinate_annihilator_left_witnesses=[[one]])
        G.certify_coordinate_bounds(s, f)
        assert G.polynomial_hom(s, f)['dimension'] == dimension
        abstract_controls.append(label)
    nf = G.Field(dict(kind='number_field', minimal_polynomial='a^3-2'))
    operator = [[[1, 1], '1'], [[0, 0], '-a']]
    s = dict(rank=1, relation_columns=1, relations=[[operator]], coordinate_annihilators=[operator],
             coordinate_annihilator_left_witnesses=[[one]])
    G.certify_coordinate_bounds(s, nf)
    assert s['polynomial_degree_bounds'] == [-1]
    abstract_controls.append('cubic coefficient-field indicial norm has no integer root')

    # A larger source can have a polynomial quotient but the ACTUAL origin jet
    # of exp(t) rejects it. This scalar example deliberately does not claim to
    # be a whole Gaussian presentation.
    false_scalar = dict(coordinate_annihilators=[[[[0, 2], '1'], [[0, 1], '-1']]],
                        polynomial_degree_bounds=[0], amplitude_basis=[[[[0, 0, 0, 0], '1']]])
    selected = G.selected_scalar_periods(false_scalar, 1, {(0,): f.K.one}, f)
    assert not selected['all_original_base_periods_polynomial']
    assert selected['certificates'][0]['original_jet_bound'] == 1
    assert selected['certificates'][0]['original_wick_jets'] == ['1', '1']
    abstract_controls.append('original derivative jet rejects spurious constant solution of exp(t)')

    rejects = []
    for kind in ['recurrence_witness', 'coordinate_witness', 'root_threshold', 'amplitude_sign',
                 'Hom_value', 'cutoff', 'original_IBP_primitive', 'actual_jet']:
        r = copy.deepcopy(results['isotropic_cubic']);s = r['exact_unlocalized_source']
        if kind == 'recurrence_witness':
            s['original_recurrence_left_witnesses'][0].append([[0]*6, '1'])
        elif kind == 'coordinate_witness':
            s['coordinate_annihilator_left_witnesses'][0][0].append([[0, 0], '1'])
        elif kind == 'root_threshold':
            s['integration_cutoff'] += 1
        elif kind == 'amplitude_sign':
            s['amplitude_basis'][0][0][1] = '-1'
        elif kind == 'Hom_value':
            r['normalized_generator_periods'][0][0][1] = '2'
        elif kind == 'cutoff':
            r['all_marker_cutoff']['slope'] = '0'
        elif kind == 'original_IBP_primitive':
            s['original_coordinate_IBP_primitives'][0][0].append([[0, 0, 0], '1'])
        else:
            r['actual_period_selection']['certificates'][0]['original_wick_jets'][0] = '2'
        try:
            G.verify_algebraic_parts(r)
        except (ValueError, AssertionError):
            rejects.append(kind)
        else:
            raise AssertionError('accepted tamper '+kind)
    for data in [dict(dimension=0, phase_terms=[]),
                 dict(dimension=1, phase_terms=[[[1], '0.5']]),
                 dict(dimension=1, field=dict(kind='number_field', minimal_polynomial='a^2-1'), phase_terms=[])]:
        try:
            G.normalize(data)
        except ValueError:
            rejects.append('invalid input '+str(data))
        else:
            raise AssertionError('accepted invalid input')
    if args.generate:
        temp = p.parents[2]/'tmp/gaussian_source_hom_cli_20261116'
        temp.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(p/'gaussian_source_hom_producer.py'),
                   str(root/'isotropic_cubic/input.json'), str(temp/'result.json'), '--workdir', str(temp/'cas')]
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=30)
        assert G.verify_algebraic_parts(json.loads((temp/'result.json').read_text()))
    receipt = dict(completed=True, canonical_proof_sha256=sha(proof),
                   generator_sha256=sha(p/'gaussian_source_hom_producer.py'), checker_sha256=sha(Path(__file__)),
                   runtime_sha256=sha(p/'gaussian_source_hom_runtime.json'),
                   source_controls=controls, abstract_hom_controls=abstract_controls,
                   original_marker_cutoff_checks=actual_checks, rejection_controls=rejects,
                   CLI_executed_in_initial_generation=True, retained_example_files_sha256=files,
                   scope='finite controls of a general algorithm; no universal theorem inferred from these controls',
                   assurance='CAS integration completeness is not independently certified or needed for the decision: original IBP, actual origin jets, complete scalar indicial bounds and all-marker recurrence replayed exactly')
    (p/'gaussian_source_hom_producer_verification.json').write_text(dump(receipt))
    print(json.dumps(dict(completed=True, source_controls=len(controls), abstract_controls=len(abstract_controls),
                          marker_checks=actual_checks, rejections=len(rejects), retained_files=len(files))))


if __name__ == '__main__':
    main()
