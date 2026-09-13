#!/usr/bin/env sage-python
"""Narrow admission gate for the X948 fibration-diversification experiment.

Only generic equations, section words and the complete retained A1 table enter.
No specialization or point search runs here. Existing construction is reused.
"""
import argparse
import csv
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, matrix, vector
from research_runtime.store import checkpoint

CAS = Path(__file__).resolve().parent
geo = SourceFileLoader('seed_foundry_geometry', str(CAS/'parent_foundry_geometry.sage')).load_module()


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def proposal_order(folder):
    source = read(folder/'inputs/source.json')
    g = matrix(QQ, source['generic_height_gram'])
    rows = read(folder/'inputs/a1-proposals.json')['rows']
    if len(rows) != 1266 or g.det() != 948:
        raise ArithmeticError('the complete retained X948 A1 input changed')
    for row in rows:
        w = vector(ZZ, row['trace_word'])
        if w*g*w != 8:
            raise ArithmeticError('trace is not norm eight')
    known = read(folder/'inputs/prior-fibrations.json')
    R = PolynomialRing(QQ, 't')
    keys = []
    for p in known['rows']:
        inv = geo.affine_invariant(R(p['A']), R(p['B']))
        keys.append({'family': p['family'], 'invariant': inv})
    if len({p['invariant']['key'] for p in keys}) != 8:
        raise ArithmeticError('five old plus three searched X948 A1 identities required')
    old_priorities = {414, 488, 1070}
    def cost(row):
        return tuple(row[k] for k in ('group_addition_upper_bound', 'support_count',
                    'maximum_absolute_coefficient', 'coefficient_l1')) + (tuple(row['trace_word']), row['priority'])
    pool = sorted((r for r in rows if r['priority'] not in old_priorities), key=cost)[:64]
    chosen = []
    while len(chosen) < 24:
        def separation(row):
            w = vector(ZZ, row['trace_word'])
            return min(16-2*abs(w*g*vector(ZZ, p['trace_word'])) for p in chosen) if chosen else 0
        best = min(pool, key=lambda r: (-separation(r), cost(r)))
        chosen.append(best)
        pool.remove(best)
    checkpoint(folder/'proposal-order.json', {
        'status': 'FROZEN_GEOMETRY_ONLY_ORDER', 'universe_count': 1266,
        'cheap_window': 64, 'maximum_constructions': 24, 'required_admitted': 8,
        'rows': chosen, 'prior_invariants': keys,
        'rule': 'First 64 unsearched trace words in retained arithmetic-construction cost order; '
                'greedy maximum minimum distance modulo sign in the source height lattice, '
                'ties by cost and literal coordinates. First eight passing exact admission gates.',
        'boundary': 'Trace-word diversity is a scheduling proxy, not a certificate of distinct O(NS) U-orbits.'})


def admission(path):
    p = read(path)
    geo.verify(path)
    g = matrix(QQ, p['source_generic_gram'])
    w = vector(ZZ, p['trace_word'])
    # Coordinates in the saturated source NS = U + MW17(-1).
    old_o = vector(ZZ, [-1, 1]+[0]*17)
    fibre = old_o + geo.source_class(w, g)
    zero = geo.source_class(vector(ZZ, p['zero_word']), g)
    sections = [geo.source_class(vector(ZZ, v), g) for v in p['section_words']]
    change = matrix(QQ, [fibre, zero, old_o]+sections)
    if change.denominator() != 1 or abs(change.det()) != 1:
        raise ArithmeticError('displayed sections do not generate the full integral NS lattice')
    if geo.pairing(fibre, fibre, g) or geo.pairing(zero, zero, g) != -2 or geo.pairing(fibre, zero, g) != 1:
        raise ArithmeticError('primitive marked U identities failed')
    if geo.pairing(old_o, old_o, g) != -2 or geo.pairing(old_o, fibre, g) or geo.pairing(old_o, zero, g):
        raise ArithmeticError('nonidentity I2 component differs')
    if any(geo.pairing(v, fibre, g) != 1 for v in sections):
        raise ArithmeticError('a purported section is not degree one')
    R = PolynomialRing(QQ, 't')
    A, B = R(p['raw_A']), R(p['raw_B'])
    delta = 4*A**3+27*B**2
    if [A.degree(), B.degree(), delta.degree()] != [8, 12, 22] or not delta.is_squarefree() or A.gcd(delta).degree():
        raise ArithmeticError('semistable fibre configuration not certified')
    if QQ(p['generic_height_gram_determinant']) != 474:
        raise ArithmeticError('saturated height determinant differs')
    if p['compactification']['after_bits'] > 256:
        raise ArithmeticError('frozen compact coefficient gate exceeded')
    certificate = {
        'status': 'PASS_SATURATED_MARKED_A1_MW16', 'parent_sha256': sha(path),
        'family': p['family'], 'generic_rank': 16, 'height_determinant': '474',
        'source_NS_determinant': '948', 'source_picard_rank': 19,
        'NS_basis_change': [list(map(str, v)) for v in change.rows()],
        'NS_basis_change_determinant': str(change.det()),
        'marked_U': [list(map(str, fibre)), list(map(str, zero+fibre))],
        'nonidentity_component': list(map(str, old_o)),
        'source_binding': p['source_input_sha256'],
        'saturation_argument': 'F,O,C and the sixteen sections generate the saturated source NS integrally. '
            'Quotient by F,O,C gives the full MW group. The single A1 root is primitive; '
            'the semistable K3 height formula also excludes torsion (height at least 4-1/2).',
        'boundary': 'Arithmetic and geometric generic rank exactly sixteen. '
            'This does not certify specialized independence or distinct abstract U-orbits.'}
    checkpoint(path.parent/'admission.json', certificate)
    return certificate


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['order', 'construct', 'verify'])
    parser.add_argument('--folder', type=Path, required=True)
    parser.add_argument('--index', type=int)
    args = parser.parse_args()
    folder = args.folder.resolve()
    if args.mode == 'order':
        proposal_order(folder)
    else:
        row = read(folder/'proposal-order.json')['rows'][args.index]
        dest = folder/'constructions'/f"{args.index:02d}-{row['priority']:05d}"
        if args.mode == 'construct':
            geo.ATLAS = folder/'inputs/native-atlas.json'
            geo.TABLE = folder/'inputs/norm8-table.tsv'
            geo.construct(row['priority'], dest)
        admission(dest/'parent.json')
