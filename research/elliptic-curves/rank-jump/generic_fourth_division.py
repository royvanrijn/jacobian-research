#!/usr/bin/env python3
"""Exact affine cochain controls and application of a uniform level-four theorem."""
from collections import deque
import argparse
import itertools
from pathlib import Path
import retrospective as r

PRIOR = r.OUT / 'rank_jump_four_division_class_creation_v1.json'
VERIFIED = r.OUT / 'rank_jump_four_division_class_creation_verification_v1.json'
PANEL = r.OUT / 'rank_jump_fresh_governing_panel_inputs_v1.json'
PROTOCOL = Path(__file__).with_name('GENERIC_FOURTH_DIVISION_PROTOCOL.json')
OUT = r.OUT / 'rank_jump_generic_fourth_division_v1.json'
LINEAR = [(0,1,3,0), (1,1,0,1), (3,0,0,1)]
I = (1,0,0,1)


def affine_cochains(m):
    zero = (0,) * (2*m)
    gens = [(zero, A) for A in LINEAR]
    gens += [(tuple(int(k==j) for k in range(2*m)), I) for j in range(2*m)]
    identity = (zero, I)
    expr = {identity: (0,0)}
    pending = deque([identity])
    piv = {}
    edges = 0
    while pending:
        t, A = pending.popleft()
        for j, (u, B) in enumerate(gens):
            C = tuple(sum(A[2*i+k]*B[2*k+l] for k in range(2)) % 4
                      for i in range(2) for l in range(2))
            v = tuple((t[2*h+i] + sum(A[2*i+k]*u[2*h+k] for k in range(2))) % 4
                      for h in range(m) for i in range(2))
            nxt = (v, C)
            e = tuple(expr[(t,A)][i] ^ sum((A[2*i+k]%2) << (2*j+k) for k in range(2))
                      for i in range(2))
            if nxt not in expr:
                expr[nxt] = e
                pending.append(nxt)
            else:
                for a, b in zip(expr[nxt], e):
                    c = a ^ b
                    while c:
                        k = c.bit_length()-1
                        if k not in piv:
                            piv[k] = c
                            break
                        c ^= piv[k]
            edges += 1
    assert len(expr) == 96 * 16**m
    zdim = 2*len(gens) - len(piv)
    assert zdim == m+3
    # Verify explicit independent coordinate cocycles, the retained derivative
    # cocycle on GL2(Z/4), and both coboundaries span the entire solution space.
    coordinates = [sum(1 << (2*(3+2*h+i)+i) for i in range(2)) for h in range(m)]
    cob = []
    for v in [(1,0), (0,1)]:
        cob.append(sum(((sum(A[2*i+k]*v[k] for k in range(2)) - v[i]) % 2) << (2*j+i)
                       for j, (_, A) in enumerate(gens) for i in range(2)))
    derivative = next(x for x in r.read(PRIOR)['abstract_full_group_controls'][0]['cocycle_generator_assignments']
                      if r.reduce(x, r.basis(cob)))
    full = cob + coordinates + [derivative]
    assert r.rank(full) == zdim
    assert all((a & c).bit_count()%2 == 0 for a in full for c in piv.values())
    return {'generic_points': m, 'order': len(expr), 'cayley_edges': edges,
            'constraint_rank': len(piv), 'Z1_dimension': zdim, 'B1_dimension': 2,
            'H1_dimension': zdim-2, 'coordinate_cocycles': coordinates,
            'derivative_cocycle': derivative, 'coboundaries': cob}


def independent_module_checks():
    # All subgroups of (Z/4)^2, by exhaustive adjoining of cyclic generators.
    elements = list(itertools.product(range(4), repeat=2))
    zero = frozenset([(0,0)])
    seen = {zero}
    pending = [zero]
    while pending:
        H = pending.pop()
        for v in elements:
            J = frozenset(((h[0]+a*v[0])%4, (h[1]+a*v[1])%4)
                          for h in H for a in range(4))
            if J not in seen:
                seen.add(J)
                pending.append(J)
    full_reductions = [H for H in seen if len({(x%2,y%2) for x,y in H}) == 4]
    assert len(full_reductions) == 1 and len(full_reductions[0]) == 16
    commuting = []
    for B in itertools.product(range(2), repeat=4):
        def mul(X,Y):
            return tuple(sum(X[2*i+k]*Y[2*k+j] for k in range(2))%2
                         for i in range(2) for j in range(2))
        if all(mul(A,B)==mul(B,A) for A in LINEAR): commuting.append(B)
    assert commuting == [(0,0,0,0), I]
    return {'subgroups_Z4_squared': len(seen), 'subgroups_with_full_mod2_image': 1,
            'equivariant_maps_Z4_squared_to_V': [list(B) for B in commuting],
            'Hom_dimension_per_generic_point': 1}


def calculate():
    old = r.read(PRIOR)
    verified = r.read(VERIFIED)
    assert old['status'] == verified['status'] == 'PASS'
    # Check every direct binding in the retained arithmetic certificates.
    for obj in [old, verified]:
        for path, digest in obj['bindings'].items():
            assert r.digest((r.ROOT/path).read_bytes()) == digest
    panel = {c['token']: c for c in r.read(PANEL)['cases']}
    rows = []
    for row in old['rows']:
        assert row['status'] == 'PASS' and row['four_division_Selmer_dimension'] == 0
        m = len(panel[row['token']]['generic_sections'])
        assert m in (16,17)
        rows.append({'token': row['token'], 'generic_dimension': m,
            'derivative_excluding_place': row['witness']['place'],
            'translation_group_over_Q_E4': f'(Z/4)^{2*m}',
            'relative_field_degree': str(16**m),
            'translation_image_collapse': False,
            'split_Kummer_dimension': m+1, 'split_Selmer_dimension': m,
            'new_Selmer_dimensions_from_generic_fourth_division': 0,
            'proof': 'uniform theorem; actual linear mod4 image need not be full'})
    assert len(rows) == 16
    return {'schema': 'rank-jump.generic-fourth-division.v1', 'status': 'PASS',
            'affine_controls': [affine_cochains(m) for m in (0,1,2)],
            'independent_module_checks': independent_module_checks(), 'rows': rows,
            'bindings': {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
                         for p in [Path(__file__), PROTOCOL, PRIOR, VERIFIED, PANEL]},
            'boundary': 'Fourth division only. No new point, CT, class group, or higher-torsion-image computation.'}


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('mode', choices=['capture','check'])
    args=p.parse_args()
    result=calculate()
    if args.mode=='capture': r.write_new(OUT,result)
    else: assert result==r.read(OUT)
    print(result['status'], [(x['order'],x['H1_dimension']) for x in result['affine_controls']])
